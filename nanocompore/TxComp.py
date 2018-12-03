# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, Counter, defaultdict
import warnings


# Third party
from scipy.stats import mannwhitneyu, ks_2samp, ttest_ind, chi2
import statsmodels.discrete.discrete_model as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

# Local package
from nanocompore.common import NanocomporeError

# Init randon seed
np.random.seed(42)

def txCompare(ref_pos_list, methods=None, sequence_context=0, min_coverage=20, logger=None, ref=None, sequence_context_weights="uniform"):

    if sequence_context_weights != "uniform" and sequence_context_weights != "harmonic":
        raise NanocomporeError ("Invalid sequence_context_weights (uniform or harmonic)")

    n_lowcov = 0
    tests = set()
    for pos, pos_dict in enumerate(ref_pos_list):

        # Filter out low coverage positions
        lowcov = False
        for cond_dict in pos_dict["data"].values():
            for sample_val in cond_dict.values():
                if sample_val["coverage"] < min_coverage:
                    lowcov=True
        ref_pos_list[pos]["lowCov"]=lowcov

        # Perform stat tests if not low cov
        if lowcov:
            n_lowcov+=1
        else:
            res = dict()
            data = pos_dict['data']
            for met in methods:
                if met in ["MW", "KS", "TT"] :
                    pvalues = nonparametric_test(data, method=met)
                    res["{}_intensity_pvalue".format(met)]=pvalues[0]
                    res["{}_dwell_pvalue".format(met)]=pvalues[1]
                    tests.add("{}_intensity_pvalue".format(met))
                    tests.add("{}_dwell_pvalue".format(met))
                elif met == "GMM":
                    gmm_results = gmm_test(data, verbose=True)
                    res["GMM_pvalue"] = gmm_results[0]
                    res["GMM_model"] = gmm_results ################################# optional ?
                    tests.add("GMM_pvalue")

            # Save results in main
            ref_pos_list[pos]['txComp'] = res
    logger.debug("Skipping {} positions because not present in all samples with sufficient coverage".format(n_lowcov))

    # Combine pvalue within a given sequence context
    if sequence_context > 0:
        if sequence_context_weights == "harmonic":
            # Generate weights as a symmetrical harmonic series
            logger.debug ("Calculate harmonic weighs and cross correlation matrices by tests")
            weights = harmomic_series(sequence_context)
        else:
            weights = [1]*(2*sequence_context+1)

        # Collect pvalue lists per tests
        pval_list_dict = defaultdict(list)
        for pos_dict in ref_pos_list:
            if 'txComp' in pos_dict:
                for test in tests:
                    if test in pos_dict['txComp']:
                        pval_list_dict[test].append(pos_dict['txComp'][test])
        # Compute cross correlation matrix per test
        corr_matrix_dict = OrderedDict()
        for test in tests:
            corr_matrix_dict[test] = cross_corr_matrix (pval_list_dict[test], sequence_context)

        logger.debug ("Combine adjacent position pvalues with Hou's method positon per position")
        # Iterate over each positions in previously generated result dictionnary
        for mid_pos in range (len(ref_pos_list)):
            # Perform test only if middle pos is valid
            if not ref_pos_list[mid_pos]["lowCov"]:
                pval_list_dict = defaultdict(list)
                for pos in range (mid_pos-sequence_context, mid_pos+sequence_context+1):
                    # If any the positions is missing or any of the pvalues in the context is lowCov or NaN, consider it 1
                    if pos < 0 or pos >= len(ref_pos_list) or ref_pos_list[pos]["lowCov"]:
                        for test in tests:
                            pval_list_dict[test].append(1)
                    # else just extract the corresponding pvalue
                    else:
                        for test in tests:
                            pval_list_dict[test].append(ref_pos_list[pos]["txComp"][test])
                # Combine collected pvalues add add to dict
                for test in tests:
                    test_label = "{}_context_{}".format(test, sequence_context)
                    ref_pos_list[mid_pos]['txComp'][test_label] = combine_pvalues_hou (pval_list_dict[test], weights, corr_matrix_dict[test])

    return ref_pos_list

def nonparametric_test(data, method=None):

    if method in ["mann_whitney", "MW"]:
        stat_test = mannwhitneyu
    elif method in ["kolmogorov_smirnov", "KS"]:
        stat_test = ks_2samp
    elif method in ["t_test", "TT"]:
        stat_test = ttest_ind
    else:
        raise NanocomporeError ("Invalid statistical method name (MW, KS, ttest)")

    condition_labels = tuple(data.keys())
    if len(condition_labels) != 2:
        raise NanocomporeError("The %s method only supports two conditions" % method)
    condition1_intensity = np.concatenate([ rep['intensity'] for rep in data[condition_labels[0]].values() ])
    condition2_intensity = np.concatenate([ rep['intensity'] for rep in data[condition_labels[1]].values() ])
    condition1_dwell = np.concatenate([ rep['dwell'] for rep in data[condition_labels[0]].values() ])
    condition2_dwell = np.concatenate([ rep['dwell'] for rep in data[condition_labels[1]].values() ])

    pval_intensity = stat_test(condition1_intensity, condition2_intensity)[1]
    pval_dwell = stat_test(condition1_dwell, condition2_dwell)[1]
    return (pval_intensity, pval_dwell)

def gmm_test(data, log_dwell=True, verbose=False):
    condition_labels = tuple(data.keys())
    if len(condition_labels) != 2:
        raise NanocomporeError("gmm_test only supports two conditions")

    global_intensity = np.concatenate(([v['intensity'] for v in data[condition_labels[0]].values()]+[v['intensity'] for v in data[condition_labels[1]].values()]), axis=None)
    global_dwell = np.concatenate(([v['dwell'] for v in data[condition_labels[0]].values()]+[v['dwell'] for v in data[condition_labels[1]].values()]), axis=None)

    if log_dwell:
        global_dwell = np.log10(global_dwell)

    Y = [ condition_labels[0] for v in data[condition_labels[0]].values() for _ in v['intensity'] ] + [ condition_labels[1] for v in data[condition_labels[1]].values() for _ in v['intensity'] ]
    X = StandardScaler().fit_transform([(i, d) for i,d in zip(global_intensity, global_dwell)])
    gmm_mod = GaussianMixture(n_components=2, covariance_type="full", random_state=146)
    y_pred = gmm_mod.fit_predict(X)
    # Add one to each group to avoid empty clusters
    y_pred=np.append(y_pred, [0,0,1,1])
    Y.extend([condition_labels[0], condition_labels[1], condition_labels[0], condition_labels[1]])
    S1_counts = Counter(y_pred[[i==condition_labels[0] for i in Y]])
    S2_counts = Counter(y_pred[[i==condition_labels[1] for i in Y]])
    contingency_table = np.array([[S1_counts[0],S1_counts[1]],[S2_counts[0],S2_counts[1]]], dtype="int64")
    Y = pd.get_dummies(Y)
    Y['intercept']=1
    logit = sm.Logit(y_pred,Y[['intercept',condition_labels[1]]] )
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            fitmod=logit.fit(disp=0)
            pvalue, coef = fitmod.pvalues[1], fitmod.params[1]
        except ConvergenceWarning:
            fitmod, pvalue, coef = "NC", 1, "NC"
    # Return model, real_labels (w/ intercept), GMM clusters, contingency table
    if verbose:
        return (pvalue, coef, fitmod, contingency_table, gmm_mod)
    else:
        return (pvalue, coef)

def cross_corr_matrix(pvalues_vector, context=2):
    """Calculate the cross correlation matrix of the
        pvalues for a given context.
    """
    matrix=[]
    pvalues_vector = np.array(pvalues_vector)
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
    # Return a very small number if pvalue =0
    if combined_p_value == 0:
        combined_p_value = np.finfo(np.float).min
    return combined_p_value

def harmomic_series (sequence_context):
    weights = []
    for i in range(-sequence_context, sequence_context+1):
        weights.append(1/(abs(i)+1))
    return weights
