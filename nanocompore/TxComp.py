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
def paired_test (ref_pos_dict, method="mann_whitney", sequence_context=0, min_coverage=20, pval_adjust='fdr_bh'):

    # Predefine stat test
    if method=="mann_whitney":
        stat_test = mannwhitneyu
    elif method=="kolmogorov_smirnov":
        stat_test = ks_2samp
    elif method=="t_test":
        stat_test = ttest_ind
    else:
        raise NanocomporeError ("Invalid statistical method name (MW, KS, ttest)")

    # Perform pair comparison position per position if coverage is sufficient
    pos_list = []
    pval_median_list = []
    pval_dwell_list = []
    for pos, pos_dict in ref_pos_dict.items():
        if pos_dict["S1_count"] >= min_coverage and pos_dict["S2_count"] >= min_coverage:
            pos_list.append (pos)
            pval_median_list.append (stat_test (pos_dict["S1_median"], pos_dict["S2_median"])[1])
            pval_dwell_list.append (stat_test (pos_dict["S1_dwell"], pos_dict["S2_dwell"])[1])

    # Perform multiple tests correction ####################################################################### To be done at the end after all as been calculated.
    if pval_adjust:
        pval_median_list = multipletests (np.array (pval_median_list), method=pval_adjust) [1]
        pval_dwell_list = multipletests (np.array (pval_dwell_list), method=pval_adjust) [1]

    # Write results to res_dict
    res_dict = OrderedDict ()
    for pos, pval_median, pval_dwell in zip (pos_list, pval_median_list, pval_dwell_list):
        res_dict [pos] = OrderedDict()
        res_dict [pos]["median"] = pval_median
        res_dict [pos]["dwell"] = pval_dwell

    if not sequence_context:
        return res_dict

    # If a sequence context is required combine adjacent pvalues with fishers method when possible
    else:
        res_dict_combined = OrderedDict()

        for mid_pos in res_dict.keys():
            pval_median_list = []
            pval_dwell_list = []

            try:
                for pos in range (mid_pos-sequence_context, mid_pos+sequence_context+1):
                    pval_median_list.append (res_dict [pos]["median"])
                    pval_dwell_list.append (res_dict [pos]["dwell"])

                res_dict_combined [mid_pos] = OrderedDict()
                res_dict_combined [mid_pos]["median"] = stats.combine_pvalues (pval_median_list, method='fisher') [1]
                res_dict_combined [mid_pos]["dwell"] = stats.combine_pvalues (pval_dwell_list, method='fisher') [1]

            # In case at least one of the adjacent position is missing
            except KeyError:
                pass

        return res_dict_combined
