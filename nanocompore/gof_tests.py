import numpy as np

from scipy.stats import kstest
from scipy.stats import multivariate_normal
from sklearn.preprocessing import StandardScaler


def gof_test_singlerep(data, conditions, config):
    control_label = config.get_depleted_condition()
    control_label_id = config.get_condition_ids()[control_label]

    control_data = data[conditions == control_label_id, :]
    test_data = data[conditions != control_label_id, :]

    scaler = StandardScaler()
    x_control = scaler.fit_transform(control_data)
    x_control = x_control[~np.isnan(x_control)]
    x_test = scaler.transform(test_data)
    x_test = x_test[~np.isnan(x_test)]

    mean_control = np.mean(x_control, axis=0)
    cov_control = np.cov(x_control, rowvar=0)

    probs_control = multivariate_normal.pdf(x_control,
                                            mean=mean_control,
                                            cov=cov_control,
                                            allow_singular=True)
    probs_test = multivariate_normal.pdf(x_test,
                                         mean=mean_control,
                                         cov=cov_control,
                                         allow_singular=True)

    result = kstest(probs_test, probs_control)
    return result.pvalue


def gof_test_multirep(data, samples, conditions, config):
    control_label = config.get_depleted_condition()
    test_label = config.get_test_condition()

    internal_likelihoods = []
    external_likelihoods = []

    control_samples = config.get_condition_samples(control_label)
    test_samples = config.get_condition_samples(test_label)

    for sample in control_samples:
        sample_data = data[samples == sample]
        sample_mean = np.mean(sample_data, axis=0)
        sample_cov = np.cov(sample_data, rowvar=False)

        tech_replicates = [s for s in control_samples if s != sample]
        for replicate in tech_replicates:
            replicate_data = data[samples == replicate]
            likelihood = multivariate_normal.pdf(replicate_data,
                                                 mean=sample_mean,
                                                 cov=sample_cov,
                                                 allow_singular=True)
            internal_likelihoods.append(likelihood)

        for bio_replicate in test_samples:
            replicate_data = data[samples == bio_replicate]
            likelihood = multivariate_normal.pdf(replicate_data,
                                                 mean=sample_mean,
                                                 cov=sample_cov,
                                                 allow_singular=True)
            external_likelihoods.append(likelihood)

    internal_likelihoods = np.concatenate(internal_likelihoods)
    external_likelihoods = np.concatenate(external_likelihoods)

    result = kstest(external_likelihoods, internal_likelihoods)
    return result.pvalue

