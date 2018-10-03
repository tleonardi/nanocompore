# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict

# Third party
from scipy.stats import mannwhitneyu, ks_2samp, ttest_ind, combine_pvalues
from statsmodels.stats.multitest import multipletests
import numpy as np

# Local package
from nanocompore.common import NanocomporeError

#~~~~~~~~~~~~~~NON PARAMETRIC STATS METHOD~~~~~~~~~~~~~~#
def paired_test (ref_pos_dict, method="mann_whitney", sequence_context=0, min_coverage=20):

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

        ############################################################## Dirty hack to correct for high coverage bias
        # Number of batch tests dependent on coverage
        n_batch = (max((len(pos_dict["S1_median"]), len(pos_dict["S2_median"])))*2)//min_coverage
        if n_batch > 100:
            n_batch = 100

        # Compute pvalues
        for var in ("median", "dwell"):

            s1_data = pos_dict["S1_"+var]
            s2_data = pos_dict["S2_"+var]
            pval_array = np.empty (shape=n_batch, dtype=np.float64)
            for i in range (n_batch):
                pval_array[i] = stat_test (np.random.choice (s1_data, min_coverage), np.random.choice (s2_data, min_coverage))[1]
            pval = np.median (pval_array)

            ref_pos_dict[pos]["pvalue_"+var] = pval

    # If a sequence context is required combine adjacent pvalues with fishers method when possible
    if sequence_context:
        lab_median = "pvalue_median_context={}".format(sequence_context)
        lab_dwell = "pvalue_dwell_context={}".format(sequence_context)
        for mid_pos in ref_pos_dict.keys():
            pval_median_list = []
            pval_dwell_list = []
            try:
                for pos in range (mid_pos-sequence_context, mid_pos+sequence_context+1):
                    pval_median_list.append (ref_pos_dict[pos]["pvalue_median"])
                    pval_dwell_list.append (ref_pos_dict[pos]["pvalue_dwell"])
                ref_pos_dict[mid_pos][lab_median] = combine_pvalues (pval_median_list, method='fisher') [1]
                ref_pos_dict[mid_pos][lab_dwell] = combine_pvalues (pval_dwell_list, method='fisher') [1]

            # In case at least one of the adjacent position is missing
            except KeyError:
                pass

    # Return final res
    return ref_pos_dict
