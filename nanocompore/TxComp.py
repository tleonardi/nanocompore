# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, Counter, defaultdict
import warnings

# Third party
from loguru import logger
from scipy.stats import mannwhitneyu, ttest_ind, chi2, f_oneway
from scipy.stats.mstats import ks_twosamp
import statsmodels.discrete.discrete_model as dm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

# Local package
from nanocompore.common import *


class TxComp(object):
    """Compare transcript data from two samples using statistical methods"""

    def __init__(self,
                 random_state,
                 db_samples:dict,
                 univariate_test:str,
                 fit_gmm:bool,
                 gmm_test:str,
                 sequence_context:int=0,
                 sequence_context_weighting:str="uniform", # or: "harmonic"
                 min_coverage:int=20,
                 allow_anova_warnings:bool=False):
        self._random_state = random_state

        if len(db_samples) != 2:
            raise NanocomporeError(f"Expected two experimental conditions, found {len(db_samples)}: {', '.join(db_samples.keys())}")
        # Store sample/condition broken down for processing below:
        self._cond1, self._cond2 = tuple(db_samples.keys())
        self._cond1_samples = [n for n, _ in db_samples[self._cond1]]
        self._cond2_samples = [n for n, _ in db_samples[self._cond2]]

        self._univariate_test = univariate_test
        self._fit_gmm = fit_gmm
        self._gmm_test = gmm_test
        self._min_coverage = min_coverage
        self._sequence_context = sequence_context
        if sequence_context > 0:
            if sequence_context_weighting == "harmonic":
                # Generate weights as a symmetrical harmonic series
                self._sequence_context_weights = self.__harmonic_series()
            elif sequence_context_weighting == "uniform":
                self._sequence_context_weights = [1] * (2 * self._sequence_context + 1)
            else:
                raise NanocomporeError("Invalid sequence context weighting ('uniform' or 'harmonic')")
        self._allow_anova_warnings = allow_anova_warnings

        # Pre-select univariate test:
        if self._univariate_test == "MW":
            self._stat_test = lambda x, y: mannwhitneyu(x, y, alternative='two-sided')
        elif self._univariate_test == "KS":
            self._stat_test = ks_twosamp
        elif self._univariate_test == "ST":
            self._stat_test = lambda x, y: ttest_ind(x, y, equal_var=False)
        elif not self._univariate_test:
            self._stat_test = None
        else:
            raise NanocomporeError("Invalid univariate test name (MW, KS, ST)")


    def __call__(self, ref_id, kmer_data):
        """Perform comparisons for one transcript ('ref_id') given k-mer data"""
        logger.debug("TxComp()")
        n_lowcov = 0
        results = {}
        n_univariate_tests = 0
        n_gmm_tests = 0
        for pos, pos_dict in kmer_data.items():
            logger.trace(f"Processing position {pos}")
            # Filter out low coverage positions
            if self.__has_low_coverage(pos_dict):
                logger.trace(f"Position {pos} has low coverage, skipping")
                n_lowcov += 1
                continue

            intensities = {k: v["intensity"] for k, v in pos_dict.items()}
            dwell_times = {k: v["dwelltime"] for k, v in pos_dict.items()}

            cond1_intensity = np.concatenate([intensities[rep] for rep in self._cond1_samples])
            cond2_intensity = np.concatenate([intensities[rep] for rep in self._cond2_samples])
            cond1_dwell = np.concatenate([dwell_times[rep] for rep in self._cond1_samples])
            cond2_dwell = np.concatenate([dwell_times[rep] for rep in self._cond2_samples])

            # Perform stat tests
            res = {}
            if self._univariate_test:
                logger.trace(f"Running {self._univariate_test} test on position {pos}")
                res["intensity_pvalue"] = self.__univariate_test(cond1_intensity, cond2_intensity)
                res["dwelltime_pvalue"] = self.__univariate_test(cond1_dwell, cond2_dwell)
                n_univariate_tests += 2

            if self._fit_gmm:
                logger.trace(f"Fitting GMM on position {pos}")
                try:
                    gmm_results = self.__gmm_fit(intensities, dwell_times)
                except:
                    raise NanocomporeError(f"Error running GMM test on transcript {ref_id}")
                for key, value in gmm_results.items():
                    res["gmm_" + key] = value
                if gmm_results["pvalue"] is not None: # TODO: can this be 'nan'?
                    n_gmm_tests += 1

            # Calculate shift statistics
            logger.trace(f"Calculatign shift stats for {pos}")
            res["shift_stats"] = self.__shift_stats(cond1_intensity, cond2_intensity,
                                                    cond1_dwell, cond2_dwell)
            # Save results in main
            logger.trace(f"Saving test results for {pos}")
            results[pos] = res

        logger.debug(f"Skipped {n_lowcov} positions because not present in all samples with sufficient coverage")

        if self._sequence_context > 0:
            if self._univariate_test:
                self.__combine_adjacent_pvalues(results, "intensity_pvalue")
                self.__combine_adjacent_pvalues(results, "dwelltime_pvalue")
            if self._fit_gmm and self._gmm_test:
                self.__combine_adjacent_pvalues(results, "gmm_pvalue")

        return (results, n_univariate_tests, n_gmm_tests)


    def __combine_adjacent_pvalues(self, results, pvalue_key):
        logger.debug(f"Calculating cross correlation matrix for '{pvalue_key}'")
        # Collect pvalue list for test
        pval_list = []
        for res_dict in results.values():
            # TODO: avoid 'None'/'np.nan' checks below by checking and replacing here?
            pval_list.append(res_dict.get(pvalue_key))
        # Compute cross correlation matrix
        corr_matrix = self.__cross_corr_matrix(pval_list)

        logger.debug("Combine adjacent position pvalues with Hou's method position by position")
        combined_label = f"{pvalue_key}_context"
        # Iterate over each position in previously generated result dictionary
        for mid_pos, res_dict in results.items():
            mid_pvalue = res_dict[pvalue_key]
            # If the mid p-value is missing or NaN, also set the context p-value to missing/NaN
            if (mid_pvalue is None) or np.isnan(mid_pvalue):
                results[mid_pos][combined_label] = mid_pvalue
                continue
            # Otherwise collect adjacent p-values and combine them:
            pval_list = []
            for pos in range(mid_pos - self._sequence_context, mid_pos + self._sequence_context + 1):
                # If any of the positions is missing or any of the p-values in the context is NaN, consider it 1
                if (pos not in results) or (results[pos][pvalue_key] is None) or np.isnan(results[pos][pvalue_key]):
                    pval_list.append(1)
                else: # just extract the corresponding pvalue
                    pval_list.append(results[pos][pvalue_key])
            # Combine collected pvalues and add to dict
            results[mid_pos][combined_label] = self.__combine_pvalues_hou(pval_list, corr_matrix)


    def __univariate_test(self, cond1_values, cond2_values):
        pvalue = self._stat_test(cond1_values, cond2_values)[1]
        if pvalue == 0:
            return np.finfo(np.float).tiny
        return pvalue


    def __gmm_fit(self, intensities, dwell_times):
        # Dictionary Sample_ID:Condition_label
        sample_condition_labels = dict([(k, self._cond1) for k in self._cond1_samples] +
                                       [(k, self._cond2) for k in self._cond2_samples])

        # Merge the intensities and dwell times of all samples in a single array
        global_intensity = np.concatenate(list(intensities.values()))
        global_dwell = np.concatenate(list(dwell_times.values()))

        # Generate the intensity and dwell time array
        X = np.array([(i, d) for i, d in zip(global_intensity, global_dwell)])

        # Generate an array of sample IDs
        Y = np.concatenate((np.repeat(self._cond1_samples, [len(intensities[k]) for k in self._cond1_samples]),
                            np.repeat(self._cond2_samples, [len(intensities[k]) for k in self._cond2_samples])))

        gmm_mod, gmm_type, gmm_ncomponents = self.__fit_best_gmm(X, max_components=2, cv_types=['full'])

        if gmm_ncomponents == 2:
            # Assign data points to the clusters
            y_pred = gmm_mod.predict(X)
            counters = dict()
            # Count how many reads in each cluster for each sample
            for lab in self._cond1_samples + self._cond2_samples:
                counters[lab] = Counter(y_pred[[i == lab for i in Y]])
            cluster_counts = self.__count_reads_in_cluster(counters)
            if self._gmm_test == "anova":
                pvalue, stat, details = self.__gmm_anova_test(counters, sample_condition_labels, gmm_ncomponents)
            elif self._gmm_test == "logit":
                pvalue, stat, details = self.__gmm_logit_test(Y, y_pred, sample_condition_labels)
        else:
            pvalue = stat = details = cluster_counts = None

        return {"model": gmm_mod, "cluster_counts": cluster_counts, "pvalue": pvalue, "test_stat": stat,
                "test_details": details}


    def __fit_best_gmm(self, X, max_components=2, cv_types=['spherical', 'tied', 'diag', 'full']):
       # Loop over multiple cv_types and n_components and for each fit a GMM
        # calculate the BIC and retain the lowest
        lowest_bic = np.infty
        bic = []
        n_components_range = range(1, max_components + 1)
        for cv_type in cv_types:
            for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
                gmm = GaussianMixture(n_components=n_components, covariance_type=cv_type,
                                      random_state=self._random_state)
                gmm.fit(X)
                bic.append(gmm.bic(X))
                if bic[-1] < lowest_bic:
                    lowest_bic = bic[-1]
                    best_gmm = gmm
                    best_gmm_type = cv_type
                    best_gmm_ncomponents = n_components
        return (best_gmm, best_gmm_type, best_gmm_ncomponents)


    def __gmm_anova_test(self, counters, sample_condition_labels, gmm_ncomponents):
        labels = []
        logr = []
        for sample, counter in counters.items():
            # Save the condition label the corresponds to the current sample
            labels.append(sample_condition_labels[sample])
            # The Counter dictionaries in counters are not ordered
            # The following line enforces the order and adds 1 to avoid empty clusters
            ordered_counter = [counter[i] + 1 for i in range(gmm_ncomponents)]
            total = sum(ordered_counter)
            normalised_ordered_counter = [i / total for i in ordered_counter]
            # Loop through ordered_counter and divide each value by the first
            logr.append(np.log(normalised_ordered_counter[0] / (1 - normalised_ordered_counter[0])))
        logr = np.around(np.array(logr), decimals=9)
        logr_s1 = [logr[i] for i, l in enumerate(labels) if l == self._cond1]
        logr_s2 = [logr[i] for i, l in enumerate(labels) if l == self._cond2]
        # If the SS for either array is 0, skip the anova test
        if sum_of_squares(logr_s1 - np.mean(logr_s1)) == 0 and sum_of_squares(logr_s2 - np.mean(logr_s2)) == 0:
            if not self._allow_anova_warnings:
                raise NanocomporeError("While doing the Anova test we found a sample with within variance = 0. Use --allow_anova_warnings to ignore.")
            else:
                aov_table = "Within variance is 0"
                aov_pvalue = np.finfo(np.float).tiny
        else:
            with warnings.catch_warnings():
                # Convert warnings to errors in order to catch them
                warnings.filterwarnings('error')
                try:
                    aov_table = f_oneway(logr_s1, logr_s2)
                    aov_pvalue = aov_table.pvalue
                except RuntimeWarning:
                    if not self._allow_anova_warnings:
                        raise NanocomporeError("While doing the Anova test a runtime warning was raised. Use --allow_anova_warnings to ignore.")
                    else:
                        warnings.filterwarnings('default')
                        aov_table = f_oneway(logr_s1, logr_s2)
                        aov_pvalue = np.finfo(np.float).tiny
        if aov_pvalue == 0:
            raise NanocomporeError("The Anova test returned a p-value of 0. This is most likely an error somewhere")
        # Calculate the delta log odds ratio, i.e. the difference of the means of the log odds ratios between the two conditions
        aov_delta_logit = float(np.mean(logr_s1) - np.mean(logr_s2))
        aov_details = {'table': aov_table, 'log_ratios': logr}
        return (aov_pvalue, aov_delta_logit, aov_details)


    def __gmm_logit_test(self, Y, y_pred, sample_condition_labels):
        Y = [sample_condition_labels[i] for i in Y]
        y_pred = np.append(y_pred, [0, 0, 1, 1])
        Y.extend([self._cond1, self._cond2, self._cond1, self._cond2])
        Y = pd.get_dummies(Y)
        Y['intercept'] = 1
        logit = dm.Logit(y_pred, Y[['intercept', self._cond2]])
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                logit_mod = logit.fit(disp=0)
                logit_pvalue, logit_coef = logit_mod.pvalues[1], logit_mod.params[1]
            except ConvergenceWarning:
                logit_mod, logit_pvalue, logit_coef = None, 1, None
        if logit_pvalue == 0:
            logit_pvalue = np.finfo(np.float).tiny
        logit_details = {'model': logit_mod}
        return (logit_pvalue, logit_coef, logit_details)


    @staticmethod
    def __count_reads_in_cluster(counters):
        cluster_counts = list()
        for k, v in counters.items():
            cluster_counts.append("%s:%s/%s" % (k, v[0], v[1]))
        cluster_counts = "_".join(cluster_counts)
        return cluster_counts


    @staticmethod
    def __shift_stats(condition1_intensity, condition2_intensity, condition1_dwell, condition2_dwell):
        """Calculate shift statistics"""
        shift_stats = OrderedDict([
            ('c1_mean_intensity', np.mean(condition1_intensity)),
            ('c2_mean_intensity', np.mean(condition2_intensity)),
            ('c1_median_intensity', np.median(condition1_intensity)),
            ('c2_median_intensity', np.median(condition2_intensity)),
            ('c1_sd_intensity', np.std(condition1_intensity)),
            ('c2_sd_intensity', np.std(condition2_intensity)),
            ('c1_mean_dwelltime', np.mean(condition1_dwell)),
            ('c2_mean_dwelltime', np.mean(condition2_dwell)),
            ('c1_median_dwelltime', np.median(condition1_dwell)),
            ('c2_median_dwelltime', np.median(condition2_dwell)),
            ('c1_sd_dwelltime', np.std(condition1_dwell)),
            ('c2_sd_dwelltime', np.std(condition2_dwell))
        ])
        return shift_stats


    def __cross_corr_matrix(self, pvalues_vector):
        """Calculate the cross correlation matrix of the pvalues for a given context."""
        context = self._sequence_context
        if len(pvalues_vector) < (context * 3) + 3:
            raise NanocomporeError(f"Not enough p-values for a context of {context}")

        pvalues_vector = np.array([i if (i is not None) and not np.isnan(i) else 1 for i in pvalues_vector])
        if any(pvalues_vector == 0) or any(np.isinf(pvalues_vector)) or any(pvalues_vector > 1):
            raise NanocomporeError("At least one p-value is invalid")

        matrix = []
        s = pvalues_vector.size
        if all(p == 1 for p in pvalues_vector):
            return np.ones((context * 2 + 1, context * 2 + 1))

        for i in range(-context, context + 1):
            row = []
            for j in range(-context, context + 1):
                row.append(np.corrcoef((np.roll(pvalues_vector, i)[context:s - context]),
                                       (np.roll(pvalues_vector, j)[context:s - context]))[0][1])
            matrix.append(row)
        return np.array(matrix)


    def __combine_pvalues_hou(self, pvalues, cor_mat):
        """ Hou's method for the approximation for the distribution of the weighted
            combination of non-independent or independent probabilities.
            If any pvalue is nan, returns nan.
            https://doi.org/10.1016/j.spl.2004.11.028
            pvalues: list of pvalues to be combined
            cor_mat: a matrix containing the correlation coefficients between pvalues
            Test: when weights are equal and cor=0, hou is the same as Fisher
            print(combine_pvalues([0.1,0.02,0.1,0.02,0.3], method='fisher')[1])
            print(hou([0.1,0.02,0.1,0.02,0.3], [1,1,1,1,1], np.zeros((5,5))))
        """
        weights = self._sequence_context_weights
        # TODO: are the following sanity checks necessary/useful?
        if len(pvalues) != len(weights):
            raise NanocomporeError("Can't combine pvalues if pvalues and weights are not the same length.")
        if cor_mat.shape[0] != cor_mat.shape[1] or cor_mat.shape[0] != len(pvalues):
            raise NanocomporeError("The correlation matrix needs to be square, with each dimension equal to the length of the pvalued vector.")
        if all(p == 1 for p in pvalues):
            return 1
        if any((p == 0 or np.isinf(p) or p > 1) for p in pvalues):
            raise NanocomporeError("At least one p-value is invalid")

        # Covariance estimation as in Kost and McDermott (eq:8)
        # https://doi.org/10.1016/S0167-7152(02)00310-3
        cov = lambda r: (3.263*r)+(0.710*r**2)+(0.027*r**3)
        k = len(pvalues)
        cov_sum = np.float64(0)
        sw_sum = np.float64(0)
        w_sum = np.float64(0)
        tau = np.float64(0)
        for i in range(k):
            for j in range(i + 1, k):
                cov_sum += weights[i] * weights[j] * cov(cor_mat[i][j])
            sw_sum += weights[i]**2
            w_sum += weights[i]
            # Calculate the weighted Fisher's combination statistic
            tau += weights[i] * (-2 * np.log(pvalues[i]))
        # Correction factor
        c = (2 * sw_sum + cov_sum) / (2 * w_sum)
        # Degrees of freedom
        f = (4 * w_sum**2) / (2 * sw_sum + cov_sum)
        # chi2.sf is the same as 1 - chi2.cdf but is more accurate
        combined_p_value = chi2.sf(tau / c, f)
        # Return a very small number if pvalue = 0
        if combined_p_value == 0:
            combined_p_value = np.finfo(np.float).tiny
        return combined_p_value


    def __harmonic_series(self):
        weights = []
        for i in range(-self._sequence_context, self._sequence_context + 1):
            weights.append(1 / (abs(i) + 1))
        return weights


    @staticmethod
    def __sum_of_squares(x):
        """
        Square each element of the input array and return the sum
        """
        x = np.atleast_1d(x)
        return np.sum(x * x)


    def __has_low_coverage(self, pos_dict):
        for sample_dict in pos_dict.values():
            if sample_dict["coverage"] < self._min_coverage:
                return True
        return False
