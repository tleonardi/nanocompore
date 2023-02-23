#!/usr/bin/env python
def sort_list(strings):
    strings.sort(key=lambda x: (not 'pvalue' in x, x))
    return strings

def sort_list_v2(strings):
    target = ['GMM_logit_pvalue', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'GMM_cov_type', 'GMM_n_clust', 'cluster_counts', 'Logit_LOR']
    target_set = set(target)
    strings = [s for s in strings if s in target_set]
    strings.sort(key=lambda x: target.index(x) if x in target else len(target))
    return strings

headers = ['cluster_counts', 'GMM_cov_type', 'GMM_logit_pvalue', 'GMM_n_clust', 'KS_dwell_pvalue', 'KS_intensity_pvalue', 'logit_LOR']
print(headers)
print(sort_list(headers))
print(sort_list_v2(headers))

headers = ['KS_dwell_pvalue', 'KS_intensity_pvalue']
print(headers)
print(sort_list(headers))
print(sort_list_v2(headers))