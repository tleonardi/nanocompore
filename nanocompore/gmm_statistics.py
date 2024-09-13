from collections import Counter
import warnings, sys

from scipy.stats import f_oneway
import statsmodels.discrete.discrete_model as dm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

from nanocompore.common import *


def gmm_test(kmer_data, transcript, random_state, anova=True, logit=False, verbose=True, allow_warnings=False):
    if len(transcript.condition_labels) != 2:
        #sys.stderr.write("gmm_test only supports two conditions\n")
        raise NanocomporeError("gmm_test only supports two conditions")

    # Scale the intensity and dwell time arrays
    X = StandardScaler().fit_transform([(i, d) for i, d in zip(kmer_data.intensity, kmer_data.dwell)])

    gmm_fit = fit_best_gmm(X, max_components=2, cv_types=['full'], random_state=random_state)
    gmm_mod, gmm_type, gmm_ncomponents = gmm_fit

    if gmm_ncomponents > 2:
        raise NanocomporeError("GMM models with n_component>2 are not supported")

    aov_results = {'pvalue': np.nan, 'delta_logit': np.nan, 'table': "NC", 'cluster_counts': "NC"}
    logit_results = {'pvalue': np.nan, 'coef': np.nan, 'model': "NC"}
    cluster_counts = "NC"

    # If the best GMM has 2 clusters do an anova test on the log odd ratios
    if gmm_ncomponents == 2:
        # Assign data points to the clusters
        y_pred = gmm_mod.predict(X)
        counters = dict()
        # Count how many reads in each cluster for each sample
        for lab in transcript.sample_labels:
            counters[lab] = Counter(y_pred[[i==lab for i in kmer_data.sample_labels]])
        cluster_counts = count_reads_in_cluster(counters)

        if anova:
            aov_results = gmm_anova_test(counters,
                                         transcript.condition_labels,
                                         gmm_ncomponents,
                                         allow_warnings)

        if logit:
            logit_results = gmm_logit_test(y_pred,
                                           kmer_data,
                                           transcript.condition_labels)

    return({'anova': aov_results,
            'logit': logit_results,
            'gmm': {'GMM_cov_type': gmm_type,
                    'n_componets': gmm_ncomponents,
                    'model': gmm_mod,
                    'cluster_counts': cluster_counts}})


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


def gmm_anova_test(counters, condition_labels, gmm_ncomponents, allow_warnings=False):
    labels= []
    logr = []
    allow_warnings = True
    for sample,counter in counters.items():
        # Save the condition label the corresponds to the current sample
        labels.append(sample)
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
            aov_pvalue = np.finfo(float).tiny
    else:
        with warnings.catch_warnings():
            # Convert warnings to errors in order to catch them
            warnings.filterwarnings('error')
            #try:
            aov_table = f_oneway(logr_s1, logr_s2)
            aov_pvalue = aov_table.pvalue
            #except RuntimeWarning:
                #if not allow_warnings:
                #    raise NanocomporeError("While doing the Anova test a runtime warning was raised. Use --allow_warnings to ignore.")
                #else:
                #    warnings.filterwarnings('default')
                #    aov_table = f_oneway(logr_s1, logr_s2)
                #    aov_pvalue = np.finfo(float).tiny
    if aov_pvalue == 0:
        sys.stderr.write("The Anova test returned a p-value of 0. This is most likely an error somewhere\n")
        #raise NanocomporeError("The Anova test returned a p-value of 0. This is most likely an error somewhere")
    # Calculate the delta log odds ratio, i.e. the difference of the means of the log odds ratios between the two conditions
    aov_delta_logit=float(np.mean(logr_s1)-np.mean(logr_s2))
    aov_results = {'pvalue': aov_pvalue, 'delta_logit': aov_delta_logit, 'table': aov_table, 'log_ratios':logr}
    return(aov_results)


def gmm_logit_test(y_pred, kmer_data, condition_labels):
    Y = kmer_data.condition_labels
    y_pred=np.append(y_pred, [0, 0, 1, 1])
    Y.extend([condition_labels[0], condition_labels[1], condition_labels[0], condition_labels[1]])
    Y = pd.get_dummies(Y)
    Y['intercept'] = 1
    X = Y[['intercept', condition_labels[0]]].astype(int)
    logit = dm.Logit(y_pred, X)
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            logit_mod=logit.fit(disp=0)
            logit_pvalue, logit_coef = logit_mod.pvalues.iloc[1], logit_mod.params.iloc[1]
        except ConvergenceWarning:
            logit_mod, logit_pvalue, logit_coef = "NC", 1, "NC"
    if logit_pvalue == 0:
        logit_pvalue = np.finfo(float).tiny
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
