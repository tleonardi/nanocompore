# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, Counter


# Third party
from scipy.stats import mannwhitneyu, ks_2samp, ttest_ind, combine_pvalues, chi2_contingency
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import numpy as np

# Local package
from nanocompore.common import NanocomporeError

#~~~~~~~~~~~~~~NON PARAMETRIC STATS METHOD~~~~~~~~~~~~~~#
def paired_test (ref_pos_dict, method="mann_whitney", sequence_context=0, min_coverage=20):
    if type(sequence_context) is not int:
        raise NanocomporeError("sequence_context is not of type int")
    np.random.seed (42)
    # Predefine stat test
    if method in ["mann_whitney", "MW"]:
        method = "mann_whitney"
        stat_test = mannwhitneyu
    elif method in ["kolmogorov_smirnov", "KS"]:
        method = "kolmogorov_smirnov"
        stat_test = ks_2samp
    elif method in ["t_test", "TT"]:
        method = "t_test"
        stat_test = ttest_ind
    else:
        raise NanocomporeError ("Invalid statistical method name (MW, KS, ttest)")

    # Perform pair comparison position per position if coverage is sufficient
    for pos, pos_dict in ref_pos_dict.items():
        # Compute pvalues
        for var in ("median", "dwell"):
            s1_data = pos_dict["S1_"+var]
            s2_data = pos_dict["S2_"+var]
            pval = stat_test (s1_data, s2_data)[1]
            ref_pos_dict[pos]["pvalue_"+method+"_"+var] = pval

    # If a sequence context is required combine adjacent pvalues with fishers method when possible
    if sequence_context:
        lab_median = "pvalue_{}_median_context={}".format(method, sequence_context)
        lab_dwell = "pvalue_{}_dwell_context={}".format(method, sequence_context)
        for mid_pos in ref_pos_dict.keys():
            pval_median_list = []
            pval_dwell_list = []
            try:
                for pos in range (mid_pos-sequence_context, mid_pos+sequence_context+1):
                    pval_median_list.append (ref_pos_dict[pos]["pvalue_"+method+"_median"])
                    pval_dwell_list.append (ref_pos_dict[pos]["pvalue_"+method+"_dwell"])
                ref_pos_dict[mid_pos][lab_median] = combine_pvalues (pval_median_list, method='fisher') [1]
                ref_pos_dict[mid_pos][lab_dwell] = combine_pvalues (pval_dwell_list, method='fisher') [1]

                if ref_pos_dict[mid_pos][lab_median] == 0:
                    ref_pos_dict[mid_pos][lab_median] = np.finfo(np.float).min
                if ref_pos_dict[mid_pos][lab_dwell] == 0:
                    ref_pos_dict[mid_pos][lab_dwell] = np.finfo(np.float).min

            # In case at least one of the adjacent position is missing
            except KeyError:
                pass

    # Return final res
    return ref_pos_dict

def kmeans_test(ref_pos_dict, method="kmeans", sequence_context=0, min_coverage=5):
    if type(sequence_context) is not int:
        raise NanocomporeError("sequence_context is not of type int")
    for pos, pos_dict in ref_pos_dict.items():
        median=np.concatenate((pos_dict['S1_median'], pos_dict['S2_median']))
        dwell=np.concatenate((pos_dict['S1_dwell'], pos_dict['S2_dwell']))
        if len(pos_dict['S1_median']) != len(pos_dict['S1_dwell']) or len(pos_dict['S2_median']) != len(pos_dict['S2_dwell']):
            raise NanocomporeError ("Median and dwell time arrays have mismatiching lengths")
        Y = ["S1" for _ in pos_dict['S1_median']] + ["S2" for _ in pos_dict['S2_median']]
        X = StandardScaler().fit_transform([(m, d) for m,d in zip(median, dwell)])
        y_pred = KMeans(n_clusters=2, random_state=146).fit_predict(X)
        S1_counts = Counter(y_pred[[i=="S1" for i in Y]])
        S2_counts = Counter(y_pred[[i=="S2" for i in Y]])
        f_obs = np.array([[S1_counts[0],S1_counts[1]],[S2_counts[0],S2_counts[1]]], dtype="int64")
        f_obs
        if any([k<min_coverage for i in f_obs for k in i ]):
            pval=1
        else:
            try:
                chi = chi2_contingency(f_obs)
                pval = chi[1]
            except:
                pval=1  
        ref_pos_dict[pos]["pvalue_kmeans"] = pval

    if sequence_context:
        lab = "pvalue_kmeans_context={}".format(sequence_context)
        for mid_pos in ref_pos_dict.keys():
            pval_list = []
            try:
                for pos in range(mid_pos-sequence_context, mid_pos+sequence_context+1):
                    pval_list.append(ref_pos_dict[pos]["pvalue_kmeans"])
                ref_pos_dict[mid_pos][lab] = combine_pvalues(pval_list, method='fisher')[1]
                if ref_pos_dict[mid_pos][lab] == 0:
                   ref_pos_dict[mid_pos][lab] = np.finfo(np.float).min
            except KeyError:
                pass

    return ref_pos_dict



