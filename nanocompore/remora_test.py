#!/usr/bin/env python3
from __future__ import print_function
import re

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
warnings.simplefilter('ignore', ConvergenceWarning)

import pod5
import pysam
import numpy as np
import math
import pandas as pd
import itertools
from tqdm import tqdm
from scipy.stats import ttest_ind

from pkg_resources import resource_filename

from remora import io, refine_signal_map, util



from collections import OrderedDict, Counter, defaultdict
import sys
# Third party
from scipy.stats.mstats import ks_twosamp
import statsmodels.discrete.discrete_model as dm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture


def gmm_test(pos, intensity, dwell, labels, sample_2_condition, random_state=26, anova=True, logit=False, verbose=True, allow_warnings=False):
    # Condition labels
    condition_labels = labels
    # List of sample labels
    sample_labels = labels

    #Forcing conditionals for testing
    anova = False
    logit = True

    # Scale the intensity and dwell time arrays
    X = StandardScaler().fit_transform([(i, d) for i,d in zip(intensity, dwell)])

    # Generate an array of sample labels
    Y = labels

    gmm_fit = fit_best_gmm(pos, X, max_components=2, cv_types=['full'], random_state=random_state)
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
            aov_results = gmm_anova_test(counters, sample_2_condition, condition_labels, gmm_ncomponents, allow_warnings)
        else:
            aov_results=None

        if logit:
            logit_results = gmm_logit_test(Y, y_pred, sample_2_condition, condition_labels)
        else:
            logit_results=None
    elif gmm_ncomponents == 1:
        aov_results = {'pvalue': np.nan, 'delta_logit': np.nan, 'table': "NC", 'cluster_counts': "NC"}
        logit_results = {'pvalue': np.nan, 'coef': np.nan, 'model': "NC"}
        cluster_counts = "NC"

    return({'anova':aov_results, 'logit': logit_results, 'gmm':{'GMM_cov_type':gmm_type, 'n_componets':gmm_ncomponents, 'model': gmm_mod, 'cluster_counts': cluster_counts}})


def fit_best_gmm(pos, X, random_state, max_components=2, cv_types=['spherical', 'tied', 'diag', 'full']):
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


def gmm_anova_test(counters, sample2condition, condition_labels, gmm_ncomponents, allow_warnings=False):
    labels= []
    logr = []
    allow_warnings = True
    for sample,counter in counters.items():
        # Save the condition label the corresponds to the current sample
        labels.append(sample2condition[sample])
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
            sys.stderr.write("While doing the Anova test we found a sample with within variance = 0. Use --allow_warnings to ignore.\n")
            #raise NanocomporeError("While doing the Anova test we found a sample with within variance = 0. Use --allow_warnings to ignore.")
        else:
            aov_table = "Within variance is 0"
            aov_pvalue = np.finfo(np.float).tiny
    else:
        aov_table = f_oneway(logr_s1, logr_s2)
        aov_pvalue = aov_table.pvalue
        #except RuntimeWarning:
            #if not allow_warnings:
            #    raise NanocomporeError("While doing the Anova test a runtime warning was raised. Use --allow_warnings to ignore.")
            #else:
            #    warnings.filterwarnings('default')
            #    aov_table = f_oneway(logr_s1, logr_s2)
            #    aov_pvalue = np.finfo(np.float).tiny
    if aov_pvalue == 0:
        sys.stderr.write("The Anova test returned a p-value of 0. This is most likely an error somewhere\n")
        #raise NanocomporeError("The Anova test returned a p-value of 0. This is most likely an error somewhere")
    # Calculate the delta log odds ratio, i.e. the difference of the means of the log odds ratios between the two conditions
    aov_delta_logit=float(np.mean(logr_s1)-np.mean(logr_s2))
    aov_results = {'pvalue': aov_pvalue, 'delta_logit': aov_delta_logit, 'table': aov_table, 'log_ratios':logr}
    return(aov_results)


def gmm_logit_test(Y, y_pred, sample2condition, condition_labels):
    Y = [ sample2condition[i] for i in Y]
    y_pred=np.append(y_pred, [0,0,1,1])
    Y.extend([condition_labels[0], condition_labels[1], condition_labels[0], condition_labels[1]])
    Y = pd.get_dummies(Y)
    Y['intercept']=1
    logit = dm.Logit(y_pred,Y[['intercept',condition_labels[1]]] )
    try:
        logit_mod=logit.fit(disp=0)
        logit_pvalue, logit_coef = logit_mod.pvalues[1], logit_mod.params[1]
    except:
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

def sum_of_squares(x):
    """
    Square each element of the input array and return the sum
    """
    x = np.atleast_1d(x)
    return np.sum(x*x)



base_dir='/hps/nobackup/birney/users/logan/nanocompore_demo_data/data'

can_pod5_path = f"{base_dir}/WT_1/pod5/output.pod5"
can_bam_path = f"{base_dir}/WT_1/basecalled_pod5.bam"
mod_pod5_path = f"{base_dir}/KD_1/pod5/output.pod5"
mod_bam_path = f"{base_dir}/KD_1/basecalled_pod5.bam"

can_pod5_fh = pod5.Reader(can_pod5_path)
can_bam_fh = pysam.AlignmentFile(can_bam_path)
mod_pod5_fh = pod5.Reader(mod_pod5_path)
mod_bam_fh = pysam.AlignmentFile(mod_bam_path)

read = next(can_pod5_fh.reads())
# Get the signal data and sample rate
sample_rate = read.run_info.sample_rate
time_per_sample = 1.0/sample_rate

level_table = resource_filename("nanocompore", "rna002_5mer_levels_v1.txt")
sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    scale_iters=0,
    do_fix_guage=True,
)

ref_reg = io.RefRegion(ctg="ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|", strand="+", start=0, end=2021)

samples_metrics, all_bam_reads = io.get_ref_reg_samples_metrics(
    ref_reg,
    [(can_pod5_fh, can_bam_fh), (mod_pod5_fh, mod_bam_fh)],
    metric="dwell_trimmean_trimsd",
    sig_map_refiner=sig_map_refiner,
    max_reads=1000,
    reverse_signal=True
)

sample2Condition ={'KD':'KD', 'WT':'WT'}
sample_labels = ['WT', 'KD']
print_counter = 0
total_pos = 0
viable_pos = 0
viable_pvalues = 0
for pos in range(ref_reg.len):
    total_pos += 1
    intensity = []
    dwell = []
    labels = []
    for i, sample in enumerate(samples_metrics):
        intensity.append(sample['trimmean'][:, pos])
        dwell.append(sample['dwell'][:, pos])
        labels.append([sample_labels[i]]*len(sample['trimmean'][:, pos]))
    intensity = np.array(list(itertools.chain(*intensity)))
    dwell = np.array(list(itertools.chain(*dwell)))
    labels = list(itertools.chain(*labels))

    data = pd.DataFrame([intensity, dwell, labels]).transpose()
    data.columns = ('intensity', 'dwell', 'labels')
    data.dropna(inplace=True)

    if (data.labels.values == 'WT').sum() > 30 and (data.labels.values == 'KD').sum() > 30:
        viable_pos += 1
        data.dwell *= time_per_sample
        data.dwell = data.dwell.apply(math.log10)

        gmm_results = gmm_test(pos, data['intensity'].to_numpy(), data['dwell'].to_numpy(), data['labels'].to_list(), sample_2_condition=sample2Condition, random_state=26)
    
        if gmm_results['logit']['pvalue'] > 0:
            viable_pvalues += 1
            if print_counter < 5:
            
                print_counter += 1
                print(pos)
                print(data)
                print(data.dtypes)
                print(gmm_results)
    
print(total_pos, viable_pos, viable_pvalues)
