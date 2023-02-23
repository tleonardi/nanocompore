# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, Counter, defaultdict
from email.policy import default
import warnings, sys


# Third party
from loguru import logger
from scipy.stats import mannwhitneyu, ttest_ind, chi2, f_oneway
from scipy.stats.mstats import ks_twosamp
import statsmodels.discrete.discrete_model as dm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

# Local package
from common import *
import TranscriptObject
import gmm_test as gmm

class TxComp():
    def __init__(self,
                 num_reference_samples,
                 num_test_samples, 
                 random_state = 26,
                 methods=None,
                 sequence_context=0,
                 min_coverage=30,
                 sequence_context_weights="uniform",
                 anova=False,
                 logit=True,
                 allow_warnings=False):
        logger.debug("TxCompare object created")

        self._random_state = random_state
        self._methods = methods
        self._sequence_context = sequence_context
        self._min_coverage = min_coverage
        self._sequence_context_weights = sequence_context_weights
        self._anova = anova
        self._logit =  logit
        self._allow_warnings = allow_warnings

        if self._sequence_context_weights != "uniform" and self._sequence_context_weights != "harmonic":
            raise NanocomporeError("Invalid sequence_context_weights (uniform or harmonic)")

        # If we have less than 2 replicates in any condition skip anova and force logit method
        if num_reference_samples < 2 or num_test_samples < 2:
            self._anova = False
            self._logit = True


    def txCompare(self, transcript):
        if transcript.enoughTxCoverage():
            logger.debug(f"TxCompare starting for {transcript.name}")
            total_results_dict = defaultdict()
            tests = set()
            valid_positions = transcript.getValidPositions()
            for pos in valid_positions:
                results_dict = {}
                results_dict['kmer_seq'] = transcript.getKmerSeq(pos)
                for method in self._methods:
                    if method.upper() in ['MW', 'KS', 'TT']:
                        try:
                            pvalues = self._nonparametric_test(transcript, pos, method=method)
                        except:
                            raise NanocomporeError(f"Error doing {method} test on {transcript.name} at pos {pos}")
                        results_dict[f"{method}_intensity_pvalue"] = pvalues[0]
                        results_dict[f"{method}_dwell_pvalue"] = pvalues[1]
                        tests.add(f"{method}_intensity_pvalue")
                        tests.add(f"{method}_dwell_pvalue")
                    
                    elif method == "GMM":
                        try:
                            gmm_results = gmm.gmm_test(transcript, pos, anova=self._anova, logit=self._logit, allow_warnings=self._allow_warnings, random_state=self._random_state)
                        except:
                            raise NanocomporeError(f"Error doing GMM test on {transcript.name}")
                        #results_dict["GMM_model"] = gmm_results['gmm']
                        results_dict["cluster_counts"] = gmm_results['gmm']['cluster_counts']
                        tests.add("cluster_counts")
                        results_dict['GMM_n_clust'] = gmm_results['gmm']['n_componets']
                        tests.add('GMM_n_clust')
                        results_dict['GMM_cov_type'] = gmm_results['gmm']['GMM_cov_type']
                        tests.add('GMM_cov_type')

                        if self._anova:
                            results_dict["GMM_anova_pvalue"] = gmm_results['anova']['pvalue']
                            results_dict["GMM_anova_model"] = gmm_results['anova']['model']
                            tests.add("GMM_anova_pvalue")
                            tests.add("GMM_anova_model")
                        if self._logit:
                            results_dict["GMM_logit_pvalue"] = gmm_results['logit']['pvalue']
                            results_dict["logit_LOR"] = gmm_results['logit']['coef']
                            tests.add("GMM_logit_pvalue")
                            tests.add("logit_LOR")
                    
                
                
                # Calculate shift statistics
                #sys.stderr.write(f"Calculatign shift stats for {pos} of {transcript.name}\n")
                logger.trace(f"Calculatign shift stats for {pos} of {transcript.name}")
                shift_stats = self._shift_stats(transcript, pos)
                for stat in shift_stats:
                    results_dict[stat] = shift_stats[stat]
                # Save results in main
                #sys.stderr.write(f"Saving test results for {pos}\n")
                logger.trace(f"Saving test results for {pos}")
                
                total_results_dict[pos] = results_dict


            n_lowcov = transcript.length - len(valid_positions)
            logger.debug(f"Skipped {n_lowcov} positions for {transcript.name} because not present in all samples with sufficient coverage")

            # Combine pvalue within a given sequence context
            if self._sequence_context > 0:
                logger.debug ("Calculate weighs and cross correlation matrices by tests")
                if self._sequence_context_weights == "harmonic":
                    # Generate weights as a symmetrical harmonic series
                    weights = self._harmomic_series(self._sequence_context)
                else:
                    weights = [1]*(2*self._sequence_context+1)

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
                        corr_matrix_dict[test] = self._cross_corr_matrix(pval_list_dict[test], self._sequence_context)

                logger.debug("Combine adjacent position pvalues with Hou's method position per position")
                # Iterate over each positions in previously generated result dictionary
                for mid_pos in range(len(total_results_dict)):
                    # Perform test only if middle pos is valid
                    if mid_pos in valid_positions:
                        pval_list_dict = defaultdict(list)
                        for pos in range(mid_pos-self._sequence_context, mid_pos+self._sequence_context+1):
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
                                test_label = "{}_context_{}".format(test, self._sequence_context)
                                # If the mid p-value is.nan, force to nan also the context p-value
                                if np.isnan(total_results_dict[mid_pos][test]):
                                    total_results_dict[mid_pos][test_label] = np.nan
                                else:
                                    total_results_dict[mid_pos][test_label] = self._combine_pvalues_hou(pval_list_dict[test], weights, corr_matrix_dict[test])

            '''
            for pos in range(1, transcript.length+1):
                if pos not in total_results_dict.keys():
                    results_dict = {}
                    for test in tests.union(set(shift_stats.keys())):
                        results_dict[test] = np.nan
                    total_results_dict[pos] = results_dict
            '''            
            return total_results_dict

        else:
            logger.debug(f"Skipping {transcript.name} because there is insuffiecent coverage in all samples")
            
            
    def _nonparametric_test(self, transcript, pos, method=None):

        if method in ["mann_whitney", "MW"]:
            stat_test = lambda x,y: mannwhitneyu(x, y, alternative='two-sided')
        elif method in ["kolmogorov_smirnov", "KS"]:
            stat_test = ks_twosamp
        elif method in ["t_test", "TT"]:
            stat_test = lambda x,y: ttest_ind(x, y, equal_var=False)
        else:
            raise NanocomporeError("Invalid statistical method name (MW, KS, TT)")

        pval_intensity = stat_test(transcript.getReferenceIntensityData(pos), transcript.getTestIntensityData(pos))[1]
        if pval_intensity == 0:
            pval_intensity = np.finfo(np.float).tiny

        pval_dwell = stat_test((transcript.getReferenceDwellData(pos)), (transcript.getTestDwellData(pos)))[1]
        if pval_dwell == 0:
            pval_dwell = np.finfo(np.float).tiny

        return(pval_intensity, pval_dwell)


    def _shift_stats(self, transcript, pos):
        """Calculate shift statistics"""
        reference_intensity = transcript.getReferenceIntensityData(pos)
        test_intensity = transcript.getTestIntensityData(pos)

        reference_dwell = np.log10(transcript.getReferenceDwellData(pos))
        test_dwell = np.log10(transcript.getTestDwellData(pos))

        shift_stats = OrderedDict([
            ('c1_mean_intensity', np.mean(reference_intensity)),
            ('c2_mean_intensity', np.mean(test_intensity)),
            ('c1_median_intensity', np.median(reference_intensity)),
            ('c2_median_intensity', np.median(test_intensity)),
            ('c1_sd_intensity', np.std(reference_intensity)),
            ('c2_sd_intensity', np.std(test_intensity)),
            ('c1_mean_dwell', np.mean(reference_dwell)),
            ('c2_mean_dwell', np.mean(test_dwell)),
            ('c1_median_dwell', np.median(reference_dwell)),
            ('c2_median_dwell', np.median(test_dwell)),
            ('c1_sd_dwell', np.std(reference_dwell)),
            ('c2_sd_dwell', np.std(test_dwell))
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
            combined_p_value = np.finfo(np.float).tiny
        return combined_p_value