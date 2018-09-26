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
        stat_test = mannwhitneyu
    elif method in ["kolmogorov_smirnov", "KS"]:
        stat_test = ks_2samp
    elif method in ["t_test", "TT"]:
        stat_test = ttest_ind
    else:
        raise NanocomporeError ("Invalid statistical method name (MW, KS, ttest)")

    # Perform pair comparison position per position if coverage is sufficient
    res_dict = OrderedDict ()

    for pos, pos_dict in ref_pos_dict.items():
        if pos_dict["S1_count"] >= min_coverage and pos_dict["S2_count"] >= min_coverage:
            res_dict[pos] = OrderedDict ()

            # Number of batch tests dependent on coverage
            n_batch = (max((pos_dict["S1_count"], pos_dict["S2_count"]))*3)//min_coverage

            # Compute
            for var in ("median", "dwell"):
                s1_data = pos_dict["S1_"+var]
                s2_data = pos_dict["S2_"+var]
                pval_array = np.empty (shape=n_batch, dtype=np.float64)
                for i in range (n_batch):
                    pval_array[i] = stat_test (np.random.choice (s1_data, min_coverage), np.random.choice (s2_data, min_coverage))[1]
                res_dict[pos][var] = np.mean (pval_array)

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
