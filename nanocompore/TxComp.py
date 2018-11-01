# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, Counter
import warnings


# Third party
from scipy.stats import mannwhitneyu, ks_2samp, ttest_ind, combine_pvalues, chi2_contingency
from statsmodels.stats.multitest import multipletests
import statsmodels.discrete.discrete_model as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd
from numpy.linalg import LinAlgError

# Local package
from nanocompore.common import NanocomporeError, cross_corr_matrix, combine_pvalues_hou

def txCompare(data, methods=None, sequence_context=0, min_coverage=0, logger=None, ref=None):

    for pos in data.keys():
        condition_labels = tuple(data[pos]['data'].keys())
        if len(condition_labels) != 2:
            logger.debug("Skipping position %s of ref %s because not present in all conditions" % (pos, ref))
            continue
        # Filter out positions with low coverage
        res = dict()
        if not all( [ rep['coverage'] > min_coverage for cond in data[pos]['data'].values() for rep in cond.values() ] ):
            res['lowCov']="yes"
            continue
        else:
            res['lowCov']="no"
        for met in methods:
            if met in ["mann_whitney", "MW", "kolmogorov_smirnov", "KS", "t_test", "TT"] :
                pvalues = nonparametric_test(data[pos]['data'], method=met)
                res[met+"intensity"]=pvalues[0]
                res[met+"dwell"]=pvalues[1]
            elif met in ["gmm", "GMM"]:
                res[met] = gmm_test(data[pos]['data'], verbose=True)
            else:
                raise NanocomporeError("The method %s is not implemented" % met)
        data[pos]['txComp'] = res
    return data

def nonparametric_test(data, method=None):
    np.random.seed(42)
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
    return pval_intensity, pval_dwell

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
        except ConvergenceWarning:
            fitmod="NC"
    # Return model, real_labels (w/ intercept), GMM clusters, contingency table
    if fitmod != "NC":
        if verbose:
            return (fitmod.pvalues[1], fitmod.params[1], fitmod, contingency_table, gmm_mod)
        else:
            return (fitmod.pvalues[1], fitmod.params[1])
    else:
        return "NC"


