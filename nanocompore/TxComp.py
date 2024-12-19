from collections import Counter
from collections import OrderedDict
from collections import defaultdict

import numpy as np
import pandas as pd

from loguru import logger
from scipy.stats import mannwhitneyu, ttest_ind, chi2
from scipy.stats.mstats import ks_twosamp

import nanocompore.gmm_statistics as gmm

from nanocompore.gof_tests import gof_test_singlerep
from nanocompore.gof_tests import gof_test_multirep
from nanocompore.common import *


class TxComp():
    def __init__(self,
                 experiment,
                 config,
                 random_state=26):
        logger.debug("TxCompare object created")

        self._config = config
        self._experiment = experiment
        self._random_state = random_state

        # If we have less than 2 replicates in any condition skip anova and force logit method
        if self._any_condition_with_fewer_than_2_samples(experiment):
            logger.info("At least one condition only has 1 sample, using logit method")
            self._anova = False
            self._logit = True


    def _any_condition_with_fewer_than_2_samples(self, experiment):
        condition_counts = defaultdict(int)
        for condition in experiment.get_condition_labels():
            condition_counts[condition] = len(experiment.condition_to_samples(condition))
        if len(condition_counts) == 2:
            if any(count < 2 for cond, count in condition_counts.items()):
                return True
            else:
                return False
        else:
            raise NanocomporeError(f"There are not exactly two condition labels")


    def txCompare(self, kmer_data_list, transcript):
        logger.debug(f"TxCompare starting for {transcript.name}")
        total_results_dict = {}
        tests = set()
        valid_positions = []
        condition_label_1, condition_label_2 = self._experiment.get_condition_labels()
        for kmer_data in kmer_data_list:
            # Make sure we have sufficient reads for all conditions.
            condition_counts = Counter(kmer_data.condition_labels)
            if not all(condition_counts.get(cond, 0) >= self._config.get_min_coverage()
                       for cond in self._experiment.get_condition_labels()):
                logger.trace(f'Skipping position {kmer_data.pos} of transcript {transcript.name} due to insuffient coverage in both conditions')
                continue

            # If we have too many reads, downsample them in a way that
            # will keep the number of reads for the two conditions equal.
            max_reads = self._config.get_downsample_high_coverage()
            kmer_data = kmer_data.subsample_reads(max_reads, random_state=self._random_state)

            valid_positions.append(kmer_data.pos)
            results_dict = {}
            results_dict['kmer_seq'] = kmer_data.kmer

            comp_methods = set(self._config.get_comparison_methods())
            auto_test = None
            if 'auto' in comp_methods:
                auto_test = self._resolve_auto_test(condition_counts.values())
                tests.add("auto")
                comp_methods = comp_methods.difference({'auto'}).union({auto_test})

            for method in comp_methods:
                if method.upper() in ['MW', 'KS', 'TT']:
                    try:
                        intensity_data_1 = kmer_data.get_condition_kmer_intensity_data(condition_label_1)
                        intensity_data_2 = kmer_data.get_condition_kmer_intensity_data(condition_label_2)
                        intensity_pvalue = self._nonparametric_test(intensity_data_1, intensity_data_2, method=method)
                        results_dict[f"{method}_intensity_pvalue"] = intensity_pvalue
                        tests.add(f"{method}_intensity_pvalue")
                        if method == auto_test:
                            results_dict['auto_pvalue'] = intensity_pvalue
                            results_dict['auto_test'] = method
                    except:
                        raise NanocomporeError(f"Error doing {method} intensity test on {transcript.name} at pos {kmer_data.pos}")

                    try:
                        dwell_data_1 = kmer_data.get_condition_kmer_dwell_data(condition_label_1)
                        dwell_data_2 = kmer_data.get_condition_kmer_dwell_data(condition_label_2)
                        dwell_pvalue = self._nonparametric_test(dwell_data_1, dwell_data_2, method=method)
                        results_dict[f"{method}_dwell_pvalue"] = dwell_pvalue
                        tests.add(f"{method}_dwell_pvalue")
                    except:
                        raise NanocomporeError(f"Error doing {method} dwell test on {transcript.name} at pos {kmer_data.pos}")

                elif method == "GMM":
                    try:
                        motor_kmer = self._get_motor_kmer(kmer_data, kmer_data_list)
                        gmm_results = gmm.gmm_test(kmer_data=kmer_data,
                                                   motor_kmer=motor_kmer,
                                                   experiment=self._experiment,
                                                   anova=self._config.get_anova(),
                                                   logit=self._config.get_logit(),
                                                   allow_warnings=self._config.get_allow_warnings(),
                                                   random_state=self._random_state)
                    except:
                        raise NanocomporeError(f"Error doing GMM test on {transcript.name}")
                    relevant_gmm_results, gmm_tests = self._extract_gmm_results(gmm_results)
                    results_dict.update(relevant_gmm_results)
                    tests.update(gmm_tests)
                    if method == auto_test:
                        results_dict['auto_pvalue'] = relevant_gmm_results['GMM_logit_pvalue']
                        results_dict['auto_test'] = method
                elif method == "GOF":
                    try:
                        motor_kmer = self._get_motor_kmer(kmer_data, kmer_data_list)
                        if self._experiment.is_multi_replicate():
                            gof_pval = gof_test_multirep(kmer_data, motor_kmer=motor_kmer, experiment=self._experiment)
                        else:
                            gof_pval = gof_test_singlerep(kmer_data, motor_kmer=motor_kmer, experiment=self._experiment)
                    except:
                        raise NanocomporeError(f"Error doing GOF test on {transcript.name}")
                    results_dict[f"GOF_pvalue"] = gof_pval
                    tests.add(f"GOF_pvalue")
                    if method == auto_test:
                        results_dict['auto_pvalue'] = gof_pval
                        results_dict['auto_test'] = method

            # Calculate shift statistics
            #sys.stderr.write(f"Calculatign shift stats for {pos} of {transcript.name}\n")
            logger.trace(f"Calculatign shift stats for {kmer_data.pos} of {transcript.name}")
            c1_intensity = kmer_data.get_condition_kmer_intensity_data(condition_label_1)
            c2_intensity = kmer_data.get_condition_kmer_intensity_data(condition_label_2)
            c1_dwell = kmer_data.get_condition_kmer_dwell_data(condition_label_1)
            c2_dwell = kmer_data.get_condition_kmer_dwell_data(condition_label_2)
            shift_stats = self._shift_stats(c1_intensity=c1_intensity, c1_dwell=c1_dwell,
                                            c2_intensity=c2_intensity, c2_dwell=c2_dwell)
            for stat in shift_stats:
                results_dict[stat] = shift_stats[stat]
            # Save results in main
            #sys.stderr.write(f"Saving test results for {pos}\n")
            logger.trace(f"Saving test results for {kmer_data.pos}")

            total_results_dict[kmer_data.pos] = results_dict

        n_lowcov = transcript.length - len(valid_positions)
        logger.debug(f"Skipped {n_lowcov} positions for {transcript.name} because not present in all samples with sufficient coverage")

        # Combine pvalue within a given sequence context
        if self._config.get_sequence_context() > 0:
            logger.debug ("Calculate weighs and cross correlation matrices by tests")
            if self._config.get_sequence_context_weights() == "harmonic":
                # Generate weights as a symmetrical harmonic series
                weights = self._harmomic_series(self._config.get_sequence_context())
            else:
                weights = [1]*(2*self._config.get_sequence_context()+1)

            # Collect pvalue lists per tests
            pval_list_dict = defaultdict(list)
            for pos in range(1, transcript.length+1):
                if pos in total_results_dict.keys():
                    pos_dict = total_results_dict[pos]
                    for test in tests:
                        pval_list_dict[test].append(pos_dict[test])
                else:
                    for test in tests:
                        pval_list_dict[test].append(np.nan)


            # Compute cross correlation matrix per test
            corr_matrix_dict = OrderedDict()
            for test in tests:
                if 'pvalue' in test:
                    corr_matrix_dict[test] = self._cross_corr_matrix(pval_list_dict[test],
                                                                     self._config.get_sequence_context())

            logger.debug("Combine adjacent position pvalues with Hou's method position per position")
            # Iterate over each positions in previously generated result dictionary
            for mid_pos in range(len(total_results_dict)):
                # Perform test only if middle pos is valid
                if mid_pos in valid_positions:
                    pval_list_dict = defaultdict(list)
                    for pos in range(mid_pos-self._config.get_sequence_context(),
                                     mid_pos+self._config.get_sequence_context()+1):
                        for test in tests:
                            if 'pvalue' in test:
                                # If any of the positions are missing or any of the pvalues in the context is lowCov or NaN, consider it 1
                                if pos < 0 or pos >= len(total_results_dict) or pos not in valid_positions or np.isnan(total_results_dict[pos][test]):
                                    pval_list_dict[test].append(1)
                                # else just extract the corresponding pvalue
                                else:
                                    pval_list_dict[test].append(total_results_dict[pos][test])
                    # Combine collected pvalues and add to dict
                    for test in tests:
                        if 'pvalue' in test:
                            test_label = "{}_context_{}".format(test, self._config.get_sequence_context())
                            # If the mid p-value is.nan, force to nan also the context p-value
                            if np.isnan(total_results_dict[mid_pos][test]):
                                total_results_dict[mid_pos][test_label] = np.nan
                            else:
                                total_results_dict[mid_pos][test_label] = self._combine_pvalues_hou(pval_list_dict[test], weights, corr_matrix_dict[test])

        return total_results_dict


    def _resolve_auto_test(self, condition_counts):
        min_coverage = min(condition_counts)
        if min_coverage < 256:
            return 'KS'
        elif min_coverage < 1024:
            return 'GMM'
        elif min_coverage >= 1024:
            return 'GOF'
        return 'KS'


    def _get_motor_kmer(self, kmer, kmers):
        for other_kmer in kmers:
            if other_kmer.pos == kmer.pos + MOTOR_DWELL_EFFECT_OFFSET:
                return other_kmer
        return None


    def _extract_gmm_results(self, gmm_results):
        results = {}
        tests = set()

        #results["GMM_model"] = gmm_results['gmm']
        results["cluster_counts"] = gmm_results['gmm']['cluster_counts']
        tests.add("cluster_counts")
        results['GMM_n_clust'] = gmm_results['gmm']['n_componets']
        tests.add('GMM_n_clust')
        results['GMM_cov_type'] = gmm_results['gmm']['GMM_cov_type']
        tests.add('GMM_cov_type')

        if self._config.get_anova():
            results["GMM_anova_pvalue"] = gmm_results['anova']['pvalue']
            results["GMM_anova_model"] = gmm_results['anova']['model']
            tests.add("GMM_anova_pvalue")
            tests.add("GMM_anova_model")
        if self._config.get_logit():
            results["GMM_logit_pvalue"] = gmm_results['logit']['pvalue']
            results["logit_LOR"] = gmm_results['logit']['coef']
            tests.add("GMM_logit_pvalue")
            tests.add("logit_LOR")

        return results, tests


    def _nonparametric_test(self, data1, data2, method=None):
        if method in ["mann_whitney", "MW"]:
            stat_test = lambda x,y: mannwhitneyu(x, y, alternative='two-sided')
        elif method in ["kolmogorov_smirnov", "KS"]:
            stat_test = ks_twosamp
        elif method in ["t_test", "TT"]:
            stat_test = lambda x,y: ttest_ind(x, y, equal_var=False)
        else:
            raise NanocomporeError("Invalid statistical method name (MW, KS, TT)")

        pval = stat_test(data1, data2)[1]
        if pval == 0:
            pval = np.finfo(float).tiny

        return pval


    def _shift_stats(self, c1_intensity, c2_intensity, c1_dwell, c2_dwell):
        """Calculate shift statistics"""

        shift_stats = OrderedDict([
            ('c1_mean_intensity', np.mean(c1_intensity)),
            ('c2_mean_intensity', np.mean(c2_intensity)),
            ('c1_median_intensity', np.median(c1_intensity)),
            ('c2_median_intensity', np.median(c2_intensity)),
            ('c1_sd_intensity', np.std(c1_intensity)),
            ('c2_sd_intensity', np.std(c2_intensity)),
            ('c1_mean_dwell', np.mean(c1_dwell)),
            ('c2_mean_dwell', np.mean(c2_dwell)),
            ('c1_median_dwell', np.median(c1_dwell)),
            ('c2_median_dwell', np.median(c2_dwell)),
            ('c1_sd_dwell', np.std(c1_dwell)),
            ('c2_sd_dwell', np.std(c2_dwell))
        ])
        return(shift_stats)


    def _harmomic_series(self, sequence_context):
        weights = []
        for i in range(-sequence_context, sequence_context+1):
            weights.append(1/(abs(i)+1))
        return weights


    def _cross_corr_matrix(self, pvalues_vector, context=2):
        """ Calculate the cross correlation matrix of the
            pvalues for a given context.
        """
        if len(pvalues_vector)<(context*3)+3:
            raise NanocomporeError(F"Not enough p-values for a context of order {context}")

        pvalues_vector = np.array([ i if not np.isnan(i) else 1 for i in pvalues_vector ])
        if any(pvalues_vector==0) or any(np.isinf(pvalues_vector)) or any(pvalues_vector>1):
            raise NanocomporeError("At least one p-value is invalid")

        matrix=[]
        s=pvalues_vector.size
        if all(p==1 for p in pvalues_vector):
            return(np.ones((context*2+1, context*2+1)))

        for i in range(-context,context+1):
            row=[]
            for j in range(-context,context+1):
                row.append(np.corrcoef((np.roll(pvalues_vector,i)[context:s-context]), (np.roll(pvalues_vector,j)[context:s-context]))[0][1])
            matrix.append(row)
        return(np.array(matrix))


    def _combine_pvalues_hou(self, pvalues, weights, cor_mat):
        """ Hou's method for the approximation for the distribution of the weighted
            combination of non-independent or independent probabilities.
            If any pvalue is nan, returns nan.
            https://doi.org/10.1016/j.spl.2004.11.028
            pvalues: list of pvalues to be combined
            weights: the weights of the pvalues
            cor_mat: a matrix containing the correlation coefficients between pvalues
            Test: when weights are equal and cor=0, hou is the same as Fisher
            print(combine_pvalues([0.1,0.02,0.1,0.02,0.3], method='fisher')[1])
            print(hou([0.1,0.02,0.1,0.02,0.3], [1,1,1,1,1], np.zeros((5,5))))
        """
        if(len(pvalues) != len(weights)):
            raise NanocomporeError("Can't combine pvalues is pvalues and weights are not the same length.")
        if( cor_mat.shape[0] != cor_mat.shape[1] or cor_mat.shape[0] != len(pvalues)):
            raise NanocomporeError("The correlation matrix needs to be squared, with each dimension equal to the length of the pvalued vector.")
        if all(p==1 for p in pvalues):
            return 1
        if any((p==0 or np.isinf(p) or p>1) for p in pvalues):
            raise NanocomporeError("At least one p-value is invalid")

        # Covariance estimation as in Kost and McDermott (eq:8)
        # https://doi.org/10.1016/S0167-7152(02)00310-3
        cov = lambda r: (3.263*r)+(0.710*r**2)+(0.027*r**3)
        k=len(pvalues)
        cov_sum=np.float64(0)
        sw_sum=np.float64(0)
        w_sum=np.float64(0)
        tau=np.float64(0)
        for i in range(k):
            for j in range(i+1,k):
                cov_sum += weights[i]*weights[j]*cov(cor_mat[i][j])
            sw_sum += weights[i]**2
            w_sum += weights[i]
            # Calculate the weighted Fisher's combination statistic
            tau += weights[i] * (-2*np.log(pvalues[i]))
        # Correction factor
        c = (2*sw_sum+cov_sum) / (2*w_sum)
        # Degrees of freedom
        f = (4*w_sum**2) / (2*sw_sum+cov_sum)
        # chi2.sf is the same as 1-chi2.cdf but is more accurate
        combined_p_value = chi2.sf(tau/c,f)
        # Return a very small number if pvalue = 0
        if combined_p_value == 0:
            combined_p_value = np.finfo(float).tiny
        return combined_p_value

