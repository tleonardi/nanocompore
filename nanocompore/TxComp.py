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
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

# Local package
from nanocompore.common import *


def txCompare(
    ref_id,
    ref_pos_list,
    random_state,
    methods=None,
    sequence_context=0,
    min_coverage=20,
    ref=None,
    sequence_context_weights="uniform",
    anova=True,
    logit=False,
    allow_warnings=False):
    logger.debug("TxCompare")

    if sequence_context_weights != "uniform" and sequence_context_weights != "harmonic":
        raise NanocomporeError("Invalid sequence_context_weights (uniform or harmonic)")

    n_lowcov = 0
    tests = set()
    # If we have less than 2 replicates in any condition skip anova and force logit method
    if not all([ len(i)>1 for i in ref_pos_list[0]['data'].values() ]):
        anova=False
        logit=True
    for pos, pos_dict in enumerate(ref_pos_list):
        logger.trace(f"Processing position {pos}")
        # Filter out low coverage positions
        lowcov = False
        for cond_dict in pos_dict["data"].values():
            for sample_val in cond_dict.values():
                if sample_val["coverage"] < min_coverage:
                    lowcov=True
        ref_pos_list[pos]["lowCov"]=lowcov

        # Perform stat tests if not low cov
        if lowcov:
            logger.trace(f"Position {pos} is low coverage, skipping")
            n_lowcov+=1
        else:
            res = dict()
            data = pos_dict['data']
            condition_labels = tuple(data.keys())
            if len(condition_labels) != 2:
                raise NanocomporeError("The %s method only supports two conditions" % method)
            condition1_intensity = np.concatenate([ rep['intensity'] for rep in data[condition_labels[0]].values() ])
            condition2_intensity = np.concatenate([ rep['intensity'] for rep in data[condition_labels[1]].values() ])
            condition1_dwell = np.concatenate([ rep['dwell'] for rep in data[condition_labels[0]].values() ])
            condition2_dwell = np.concatenate([ rep['dwell'] for rep in data[condition_labels[1]].values() ])

            for met in methods:
                logger.trace(f"Running {met} test on position {pos}")
                if met in ["MW", "KS", "TT"] :
                    try:
                        pvalues = nonparametric_test(condition1_intensity, condition2_intensity, condition1_dwell, condition2_dwell, method=met)
                    except:
                        raise NanocomporeError("Error doing {} test on reference {}".format(met, ref_id))
                    res["{}_intensity_pvalue".format(met)]=pvalues[0]
                    res["{}_dwell_pvalue".format(met)]=pvalues[1]
                    tests.add("{}_intensity_pvalue".format(met))
                    tests.add("{}_dwell_pvalue".format(met))
                elif met == "GMM":
                    try:
                        gmm_results = gmm_test(data, anova=anova, logit=logit, allow_warnings=allow_warnings, random_state=random_state)
                    except:
                        raise NanocomporeError("Error doing GMM test on reference {}".format(ref_id))
                    res["GMM_model"] = gmm_results['gmm']
                    if anova:
                        res["GMM_anova_pvalue"] = gmm_results['anova']['pvalue']
                        res["GMM_anova_model"] = gmm_results['anova']
                        tests.add("GMM_anova_pvalue")
                    if logit:
                        res["GMM_logit_pvalue"] = gmm_results['logit']['pvalue']
                        res["GMM_logit_model"] = gmm_results['logit']
                        tests.add("GMM_logit_pvalue")

            # Calculate shift statistics
            logger.trace(f"Calculatign shift stats for {pos}")
            res['shift_stats'] = shift_stats(condition1_intensity, condition2_intensity, condition1_dwell, condition2_dwell)
            # Save results in main
            logger.trace(f"Saving test results for {pos}")
            ref_pos_list[pos]['txComp'] = res
    logger.debug("Skipped {} positions because not present in all samples with sufficient coverage".format(n_lowcov))

    # Combine pvalue within a given sequence context
    if sequence_context > 0:
        logger.debug ("Calculate weighs and cross correlation matrices by tests")
        if sequence_context_weights == "harmonic":
            # Generate weights as a symmetrical harmonic series
            weights = harmomic_series(sequence_context)
        else:
            weights = [1]*(2*sequence_context+1)

        # Collect pvalue lists per tests
        pval_list_dict = defaultdict(list)
        for pos_dict in ref_pos_list:
            if 'txComp' in pos_dict:
                for test in tests:
                    pval_list_dict[test].append(pos_dict['txComp'][test])
            elif pos_dict["lowCov"]:
                for test in tests:
                    pval_list_dict[test].append(np.nan)
        # Compute cross correlation matrix per test
        corr_matrix_dict = OrderedDict()
        for test in tests:
            corr_matrix_dict[test] = cross_corr_matrix(pval_list_dict[test], sequence_context)

        logger.debug("Combine adjacent position pvalues with Hou's method position per position")
        # Iterate over each positions in previously generated result dictionary
        for mid_pos in range(len(ref_pos_list)):
            # Perform test only if middle pos is valid
            if not ref_pos_list[mid_pos]["lowCov"]:
                pval_list_dict = defaultdict(list)
                for pos in range(mid_pos-sequence_context, mid_pos+sequence_context+1):
                    for test in tests:
                        # If any of the positions is missing or any of the pvalues in the context is lowCov or NaN, consider it 1
                        if pos < 0 or pos >= len(ref_pos_list) or ref_pos_list[pos]["lowCov"] or np.isnan(ref_pos_list[pos]["txComp"][test]):
                            pval_list_dict[test].append(1)
                        # else just extract the corresponding pvalue
                        else:
                            pval_list_dict[test].append(ref_pos_list[pos]["txComp"][test])
                # Combine collected pvalues and add to dict
                for test in tests:
                    test_label = "{}_context_{}".format(test, sequence_context)
                    # If the mid p-value is.nan, force to nan also the context p-value
                    if np.isnan(ref_pos_list[mid_pos]["txComp"][test]):
                        ref_pos_list[mid_pos]['txComp'][test_label] = np.nan
                    else:
                        ref_pos_list[mid_pos]['txComp'][test_label] = combine_pvalues_hou(pval_list_dict[test], weights, corr_matrix_dict[test])

    return ref_pos_list

def nonparametric_test(condition1_intensity, condition2_intensity, condition1_dwell, condition2_dwell, method=None):

    if method in ["mann_whitney", "MW"]:
        stat_test = lambda x,y: mannwhitneyu(x, y, alternative='two-sided')
    elif method in ["kolmogorov_smirnov", "KS"]:
        stat_test = ks_twosamp
    elif method in ["t_test", "TT"]:
        stat_test = lambda x,y: ttest_ind(x, y, equal_var=False)
    else:
        raise NanocomporeError("Invalid statistical method name (MW, KS, ttest)")

    pval_intensity = stat_test(condition1_intensity, condition2_intensity)[1]
    if pval_intensity == 0: 
        pval_intensity = np.finfo(np.float).tiny

    pval_dwell = stat_test(condition1_dwell, condition2_dwell)[1]
    if pval_dwell == 0: 
        pval_dwell = np.finfo(np.float).tiny
    return(pval_intensity, pval_dwell)


def gmm_test(data, random_state, anova=True, logit=False, verbose=True, allow_warnings=False):
    # Condition labels
    condition_labels = tuple(data.keys())
    # List of sample labels
    sample_labels = list(data[condition_labels[0]].keys()) + list(data[condition_labels[1]].keys())

    if len(sample_labels) != len(set(sample_labels)):
        raise NanocomporeError("Sample labels have to be unique and it looks like some are not.")

    # Dictionary Sample_label:Condition_label
    sample_condition_labels = { sk:k for k,v in data.items() for sk in v.keys() }
    if len(condition_labels) != 2:
        raise NanocomporeError("gmm_test only supports two conditions")

    # Merge the intensities and dwell times of all samples in a single array
    global_intensity = np.concatenate(([v['intensity'] for v in data[condition_labels[0]].values()]+[v['intensity'] for v in data[condition_labels[1]].values()]), axis=None)
    global_dwell = np.concatenate(([v['dwell'] for v in data[condition_labels[0]].values()]+[v['dwell'] for v in data[condition_labels[1]].values()]), axis=None)
    global_dwell = np.log10(global_dwell)

    # Scale the intensity and dwell time arrays
    X = StandardScaler().fit_transform([(i, d) for i,d in zip(global_intensity, global_dwell)])

    # Generate an array of sample labels
    Y = [ k for k,v in data[condition_labels[0]].items() for _ in v['intensity'] ] + [ k for k,v in data[condition_labels[1]].items() for _ in v['intensity'] ]

    gmm_fit = fit_best_gmm(X, max_components=2, cv_types=['full'], random_state=random_state)
    gmm_mod, gmm_type, gmm_ncomponents = gmm_fit

    # If the best GMM has 2 clusters do an anova test on the log odd ratios
    if gmm_ncomponents == 2:
        # Assign data points to the clusters
        y_pred = gmm_mod.predict(X)
        counters = dict()
        # Count how many reads in each cluster for each sample
        for lab in sample_labels:
            counters[lab] = Counter(y_pred[[i==lab for i in Y]])
        cluster_counts = count_reads_in_cluster(counters)
        if anova:
            aov_results = gmm_anova_test(counters, sample_condition_labels, condition_labels, gmm_ncomponents, allow_warnings)
        else:
            aov_results=None

        if logit:
            logit_results = gmm_logit_test(Y, y_pred, sample_condition_labels, condition_labels)
        else:
            logit_results=None

    elif gmm_ncomponents == 1:
        aov_results = {'pvalue': np.nan, 'delta_logit': np.nan, 'table': "NC", 'cluster_counts': "NC"}
        logit_results = {'pvalue': np.nan, 'coef': "NC", 'model': "NC"}
        cluster_counts = "NC"
    else:
        raise NanocomporeError("GMM models with n_component>2 are not supported")

    return({'anova':aov_results, 'logit': logit_results, 'gmm':{'model': gmm_mod, 'cluster_counts': cluster_counts}})

def fit_best_gmm(X, random_state, max_components=2, cv_types=['spherical', 'tied', 'diag', 'full']):
   # Loop over multiple cv_types and n_components and for each fit a GMM
    # calculate the BIC and retain the lowest
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, max_components+1)
    for cv_type in cv_types:
        for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
            gmm = GaussianMixture(n_components=n_components, covariance_type=cv_type, random_state=random_state)
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
                best_gmm_type = cv_type
                best_gmm_ncomponents = n_components
    return((best_gmm, best_gmm_type, best_gmm_ncomponents))

def gmm_anova_test(counters, sample_condition_labels, condition_labels, gmm_ncomponents, allow_warnings=False):
    labels= []
    logr = []
    for sample,counter in counters.items():
        # Save the condition label the corresponds to the current sample
        labels.append(sample_condition_labels[sample])
        # The Counter dictionaries in counters are not ordered
        # The following line enforces the order and adds 1 to avoid empty clusters
        ordered_counter = [ counter[i]+1 for i in range(gmm_ncomponents)]
        total = sum(ordered_counter)
        normalised_ordered_counter = [ i/total for i in ordered_counter ]
        # Loop through ordered_counter and divide each value by the first
        logr.append(np.log(normalised_ordered_counter[0]/(1-normalised_ordered_counter[0])))
    logr = np.around(np.array(logr), decimals=9)
    logr_s1 = [logr[i] for i,l in enumerate(labels) if l==condition_labels[0]]
    logr_s2 = [logr[i] for i,l in enumerate(labels) if l==condition_labels[1]]
    # If the SS for either array is 0, skip the anova test
    if sum_of_squares(logr_s1-np.mean(logr_s1)) == 0 and sum_of_squares(logr_s2-np.mean(logr_s2)) == 0:
        if not allow_warnings:
            raise NanocomporeError("While doing the Anova test we found a sample with within variance = 0. Use --allow_warnings to ignore.")
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
                if not allow_warnings:
                    raise NanocomporeError("While doing the Anova test a runtime warning was raised. Use --allow_warnings to ignore.")
                else:
                    warnings.filterwarnings('default')
                    aov_table = f_oneway(logr_s1, logr_s2)
                    aov_pvalue = np.finfo(np.float).tiny
    if aov_pvalue == 0:
        raise NanocomporeError("The Anova test returned a p-value of 0. This is most likely an error somewhere")
    # Calculate the delta log odds ratio, i.e. the difference of the means of the log odds ratios between the two conditions
    aov_delta_logit=float(np.mean(logr_s1)-np.mean(logr_s2))
    aov_results = {'pvalue': aov_pvalue, 'delta_logit': aov_delta_logit, 'table': aov_table, 'log_ratios':logr}
    return(aov_results)

def gmm_logit_test(Y, y_pred, sample_condition_labels, condition_labels):
    Y = [ sample_condition_labels[i] for i in Y]
    y_pred=np.append(y_pred, [0,0,1,1])
    Y.extend([condition_labels[0], condition_labels[1], condition_labels[0], condition_labels[1]])
    Y = pd.get_dummies(Y)
    Y['intercept']=1
    logit = dm.Logit(y_pred,Y[['intercept',condition_labels[1]]] )
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            logit_mod=logit.fit(disp=0)
            logit_pvalue, logit_coef = logit_mod.pvalues[1], logit_mod.params[1]
        except ConvergenceWarning:
            logit_mod, logit_pvalue, logit_coef = "NC", 1, "NC"
    if logit_pvalue == 0:
        logit_pvalue = np.finfo(np.float).tiny
    logit_results = {'pvalue': logit_pvalue, 'coef': logit_coef, 'model': logit_mod}
    return(logit_results)

def count_reads_in_cluster(counters):
    cluster_counts = list()
    for k,v in counters.items():
        cluster_counts.append("%s:%s/%s" % (k, v[0], v[1]))
    cluster_counts="__".join(cluster_counts)
    return(cluster_counts)

def shift_stats(condition1_intensity, condition2_intensity, condition1_dwell, condition2_dwell):
    """Calculate shift statistics"""
    shift_stats = OrderedDict([
        ('c1_mean_intensity', np.mean(condition1_intensity)),
        ('c2_mean_intensity', np.mean(condition2_intensity)),
        ('c1_median_intensity', np.median(condition1_intensity)),
        ('c2_median_intensity', np.median(condition2_intensity)),
        ('c1_sd_intensity', np.std(condition1_intensity)),
        ('c2_sd_intensity', np.std(condition2_intensity)),
        ('c1_mean_dwell', np.mean(condition1_dwell)),
        ('c2_mean_dwell', np.mean(condition2_dwell)),
        ('c1_median_dwell', np.median(condition1_dwell)),
        ('c2_median_dwell', np.median(condition2_dwell)),
        ('c1_sd_dwell', np.std(condition1_dwell)),
        ('c2_sd_dwell', np.std(condition2_dwell))
    ])
    return(shift_stats)


def cross_corr_matrix(pvalues_vector, context=2):
    """ Calculate the cross correlation matrix of the
        pvalues for a given context.
    """
    if len(pvalues_vector)<(context*3)+3:
        raise NanocomporeError("Not enough p-values for a context of order %s"%context)

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

def combine_pvalues_hou(pvalues, weights, cor_mat):
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

def harmomic_series(sequence_context):
    weights = []
    for i in range(-sequence_context, sequence_context+1):
        weights.append(1/(abs(i)+1))
    return weights

def sum_of_squares(x):
    """
    Square each element of the input array and return the sum
    """
    x = np.atleast_1d(x)
    return np.sum(x*x)
