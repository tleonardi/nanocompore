import gc
import time
from collections import defaultdict
from typing import Union

import numpy as np
import pandas as pd
import torch
from gmm_gpu.gmm import GMM
from jaxtyping import Float, Int
from loguru import logger
from numpy.lib.stride_tricks import sliding_window_view
from scipy.stats import mannwhitneyu, ttest_ind, chi2_contingency, chi2
from scipy.stats.mstats import ks_twosamp

from nanocompore.common import NanocomporeError
from nanocompore.gof_tests import gof_test_multirep
from nanocompore.gof_tests import gof_test_singlerep
from nanocompore.transcript import Transcript


INTENSITY = 0
DWELL = 1
MOTOR = 2


class TranscriptComparator:
    """
    Compare a set of positions.
    """
    def __init__(self, config, worker, random_seed=42, dtype=torch.float32):
        self._config = config
        self._random_seed = random_seed
        self._dtype = dtype
        self._motor_dwell_offset = self._config.get_motor_dwell_offset()
        self._worker = worker


    def compare_transcript(
        self,
        transcript: Transcript,
        data: Float[torch.Tensor, "positions reads vars"],
        samples: Int[torch.Tensor, "reads"],
        conditions: Int[torch.Tensor, "reads"],
        positions: Int[torch.Tensor, "positions"],
        device: str
    ) -> tuple[Transcript, Union[pd.DataFrame, None]]:
        """
        Compare the two conditions for the given transcript.

        Parameters
        ----------
        transcript : Transcript
            The transcript for which the two conditions
            will be compared.
        data : Float[torch.Tensor, "positions reads vars"]
            Torch tensor with shape (positions, reads, vars)
            containing the data for the comparison.
        samples : Int[torch.Tensor, "reads"]
            Array of sample ids with size *reads*.
        conditions : Int[torch.Tensor, "reads"]
            Array of condition ids with size *reads*.
        positions : Int[torch.Tensor, "positions"]
            Array of reference positions on the transcript
            with size *positions*.
        device : str
            String indicating the device that will be used
            for the comparison. Indicates if CPU or GPU
            should be used.

        Returns
        -------
        tuple[Transcript, pd.DataFrame]
            Tuple with the Transcript and the comparison
            results stored in a pandas DataFrame.
        """
        n_positions = len(positions)
        if n_positions == 0:
            return (transcript, None)

        results = pd.DataFrame({'transcript_id': transcript.id,
                                'pos': positions.cpu()})

        self._worker.log("trace", "Start shift stats")
        t = time.time()
        self.retry(lambda: self._add_shift_stats(results, data, conditions, device),
                   exception=torch.OutOfMemoryError)
        if device.startswith('cuda'):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
        self._worker.log("debug", f"Finished shift stats ({time.time() - t})")

        std = nanstd(data, 1).unsqueeze(1)
        outliers = (((data - data.nanmean(1, keepdim=True)) / std).abs() > 3).any(2)
        data[outliers] = np.nan

        # Standardize the data
        std = nanstd(data, 1)
        data = (data - data.nanmean(1).unsqueeze(1)) / std.unsqueeze(1)

        # If we get the same value for a variable in most of the reads
        # the others will be removed as outliers and we cannot calculate
        # the std to scale the data. I've noticed this only in Uncalled4
        # because sometimes it would return the -32768 for the intensity
        # for all/most reads.
        bad_stds = std.isnan().any(1) | (std == 0).any(1)
        if bad_stds.any():
            self._worker.log(
                    "warning",
                    "The standard deviation cannot be calculated for some positions on "
                    f"transcript {transcript.name}, but are required for scaling the data. "
                    f"The positions {positions[bad_stds].tolist()} will be skipped.")
            data = data[~bad_stds]
            positions = positions[~bad_stds]
            n_positions = positions.shape[0]
            results.drop(np.arange(results.shape[0])[bad_stds.cpu()],
                         inplace=True,
                         axis=0)

        auto_test_mask = None
        has_auto = 'auto' in self._config.get_comparison_methods()
        if has_auto:
            auto_test_mask = self._auto_test_mask(data)
        test_masks = self._get_test_masks(auto_test_mask, n_positions)

        for test, mask in test_masks.items():
            self._worker.log("debug", f"Start {test}")
            t = time.time()
            test_results = self._run_test(test,
                                          data[mask, :, :],
                                          samples,
                                          conditions,
                                          device=device)
            self._worker.log("debug", f"Finished {test} ({time.time() - t})")
            self._merge_results(results, test_results, test, mask, auto_test_mask, n_positions)

        context = self._config.get_sequence_context()
        if context > 0:
            tests = [col for col in results.columns if 'pvalue' in col]
            for test in tests:
                combined_pvals = self._combine_context_pvalues(results, test)

                label = f"{test}_context_{context}"
                results[label] = combined_pvals

        return transcript, results


    def _get_test_masks(self, auto_test_mask, n_positions):
        base_tests = set(self._config.get_comparison_methods()) - {'auto'}
        if auto_test_mask is not None:
            additional_tests = set(auto_test_mask) - base_tests
        else:
            additional_tests = []
        masks = {}
        for test in base_tests:
            masks[test] = np.ones(n_positions, dtype=bool)
        for test in additional_tests:
            masks[test] = auto_test_mask == test
        return masks


    def _auto_test_mask(self, data):
        cov = (~data[:, :, INTENSITY].isnan()).cpu().sum(1).numpy()
        return np.where(cov < 500, 'KS', 'GMM')


    def _auto_test_pvalue(self, test):
        if test == 'KS':
            return 'KS_intensity_pvalue'
        elif test == 'GMM':
            return 'GMM_chi2_pvalue'
        elif test == 'GOF':
            return 'GOF_pvalue'
        else:
            raise NotImplementedError("Unhandled test {test} as automatic test.")


    def _transcript_reads(self, kmers):
        """
        Get a list of kmers and returns a dict of type:
        transcript_id -> {read_id -> index}
        """
        transcript_reads = defaultdict(set)
        for kmer in kmers:
            transcript_reads[kmer.transcript_id].update(kmer.reads)
        return {ref_id: dict(zip(reads, range(len(reads))))
                for ref_id, reads in transcript_reads.items()}


    def _run_test(self, test, test_data, samples, conditions, device):
        if test.upper() in ['MW', 'KS', 'TT']:
            return self._nonparametric_test(test, test_data, conditions)
        elif test.upper() == 'GMM':
            return self._gmm_test(test_data, samples, conditions, device)
        elif test.upper() == 'GOF':
            return self._gof_test(test_data, samples, conditions)
        else:
            raise NotImplementedError(f"Test {test} is not implemented!")


    def _merge_results(self, results, test_results, test, mask, auto_test_mask, n_positions):
        has_auto = 'auto' in self._config.get_comparison_methods()

        if has_auto and 'auto_pvalue' not in results:
            results['auto_pvalue'] = np.full(n_positions, np.nan)
            results['auto_test'] = auto_test_mask

        for column, values in test_results.items():
            results.loc[mask, column] = values
            if has_auto and column == self._auto_test_pvalue(test):
                rows = auto_test_mask == test
                results.loc[rows, 'auto_pvalue'] = results.loc[rows, column]


    def _nonparametric_test(self, test, test_data, conditions):
        if test in ["mann_whitney", "MW"]:
            stat_test = lambda x, y: mannwhitneyu(x, y, alternative='two-sided')
        elif test in ["kolmogorov_smirnov", "KS"]:
            stat_test = ks_twosamp
        elif test in ["t_test", "TT"]:
            stat_test = lambda x, y: ttest_ind(x, y, equal_var=False)
        else:
            raise NanocomporeError("Invalid statistical method name (MW, KS, TT)")

        cond0_mask = conditions == 0

        intensity_pvals = []
        dwell_pvals = []
        for i in range(test_data.shape[0]):
            cond0_data = test_data[i, cond0_mask, :]
            cond1_data = test_data[i, ~cond0_mask, :]
            try:
                pval = stat_test(self._drop_nans(cond0_data[:, INTENSITY]).cpu(),
                                 self._drop_nans(cond1_data[:, INTENSITY]).cpu()).pvalue
                intensity_pvals.append(pval)
            except Exception as e:
                self._worker.log("error",
                                 f"Error in nonparametric test on the intensity: {e}")
                intensity_pvals.append(np.nan)
            try:
                pval = stat_test(self._drop_nans(cond0_data[:, DWELL]).cpu(),
                                 self._drop_nans(cond1_data[:, DWELL]).cpu()).pvalue
                dwell_pvals.append(pval)
            except Exception as e:
                self._worker.log("error",
                                 f"Error in nonparametric test on the dwell time: {e}")
                dwell_pvals.append(np.nan)
        return {f'{test}_intensity_pvalue': intensity_pvals,
                f'{test}_dwell_pvalue': dwell_pvals}


    def _drop_nans(self, values):
        return values[~values.isnan()]


    def _gof_test(self, test_data, samples, conditions):
        if self._config.is_multi_replicate():
            return self._gof_test_multirep(test_data, samples, conditions)
        else:
            return self._gof_test_singlerep(test_data, conditions)


    def _gof_test_singlerep(self, test_data, conditions):
        pvals = []
        for i in range(test_data.shape[0]):
            pval = gof_test_singlerep(test_data[i, :, :], conditions, self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gof_test_multirep(self, test_data, samples, conditions):
        pvals = []
        for i in range(test_data.shape[0]):
            pval = gof_test_multirep(test_data[i, :, :],
                                     samples,
                                     conditions,
                                     self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gmm_test(self, test_data, samples, conditions, device):
        # If we don't use the motor dwell time
        # we would get a tensor that only contains
        # intensity and dwell and we can directly
        # run the tests.
        if self._motor_dwell_offset == 0:
            return self._gmm_test_split(test_data, samples, conditions, device)

        # If we use the motor dwell time than we have
        # to consider some additional complications:
        # For each position some reads have motor dwell
        # and others don't. We split the positions into
        # two sets: one set of positions where all reads
        # have motor dwell time - those would be tested
        # with a 3D GMM - and a second set of positions
        # where some reads are missing the motor dwell -
        # those would be tested with a 2D GMM.
        dim3_data, dim2_data, split = self._split_by_ndim(test_data)
        columns = set()
        if dim3_data.shape[0] > 0:
            dim3_results = self._gmm_test_split(dim3_data,
                                                samples,
                                                conditions,
                                                device)
            columns.update(dim3_results.keys())
        if dim2_data.shape[0] > 0:
            dim2_results = self._gmm_test_split(dim2_data,
                                                samples,
                                                conditions,
                                                device)
            columns.update(dim2_results.keys())

        results = {}
        for column in columns:
            if dim3_data.shape[0] > 0:
                dim3_values = dim3_results[column]
            if dim2_data.shape[0] > 0:
                dim2_values = dim2_results[column]
            dim3_i = 0
            dim2_i = 0

            merged = []
            for s in split:
                if s == 0:
                    merged.append(dim3_values[dim3_i])
                    dim3_i += 1
                else:
                    merged.append(dim2_values[dim2_i])
                    dim2_i += 1
            results[column] = merged
        return results


    def _gmm_test_split(self, data, samples, conditions, device):
        test_data = data.to(self._dtype)
        s = time.time()
        def fit_model(components):
            gmm = GMM(n_components=components,
                      device=device,
                      random_seed=self._random_seed,
                      dtype=self._dtype)
            gmm.fit(test_data)
            return gmm
        gmm1 = self.retry(lambda: fit_model(1), exception=torch.OutOfMemoryError)
        bic1 = gmm1.bic(test_data)
        del gmm1
        gc.collect()
        gmm2 = self.retry(lambda: fit_model(2), exception=torch.OutOfMemoryError)

        self._worker.log("debug", f"GMM fitting time: {time.time() - s}")
        pred = gmm2.predict(test_data, force_cpu_result=False)
        pred = torch.where(test_data[:, :, 0].isnan(), np.nan, pred)

        bic2 = gmm2.bic(test_data)
        del gmm2
        gc.collect()
        if device.startswith('cuda'):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()

        contingencies = get_contigency_matrices(conditions, pred).to(device=device)
        counts = self._get_cluster_counts(contingencies, samples, conditions, pred)
        # We add 1 to all cells in all contingency
        # matrices to make sure we don't encounter
        # a devision by zero.
        contingencies += 1
        lors = calculate_lors(contingencies)
        pvals = np.array([chi2_contingency(cont).pvalue
                          for cont in contingencies.cpu().numpy()])
        # If a single component fits the data
        # better we ignore the results from
        # the two-component GMM.
        ignored_tests = (bic1 <= bic2).numpy()
        lors[ignored_tests] = np.nan
        pvals[ignored_tests] = np.nan

        return {'GMM_chi2_pvalue': pvals,
                'GMM_LOR': lors.cpu().numpy()} | counts


    def _split_by_ndim(self, test_data):
        valid = ~test_data.isnan()
        with_motor = (valid[:, :, INTENSITY] & valid[:, :, MOTOR]).sum(1)
        total = valid[:, :, INTENSITY].sum(1)
        use_motor = with_motor >= 0.9 * total

        return (test_data[use_motor],
                test_data[~use_motor, :, :MOTOR],
                (~use_motor).int())


    def _add_shift_stats(self, results, data, conditions, device):
        data = data.to(device)

        stats = {0: {}, 1: {}}
        cond0_mask = (conditions == 0)[None, :, None].to(device)

        cond0_data = torch.where(cond0_mask, data, np.nan)
        stats[0]['mean'] = cond0_data.nanmean(1).cpu()
        stats[0]['median'] = cond0_data.nanmedian(1).values.cpu()
        stats[0]['sd'] = nanstd(cond0_data, 1).cpu()

        cond1_data = torch.where(~cond0_mask, data, np.nan)
        stats[1]['mean'] = cond1_data.nanmean(1).cpu()
        stats[1]['median'] = cond1_data.nanmedian(1).values.cpu()
        stats[1]['sd'] = nanstd(cond1_data, 1).cpu()

        dims = {
                'intensity': INTENSITY,
                'dwell': DWELL,
               }
        if self._motor_dwell_offset > 0:
            dims['motor_dwell'] = MOTOR

        for cond in [0, 1]:
            for stat in ['mean', 'median', 'sd']:
                for dim, dim_index in dims.items():
                    label = f"c{cond+1}_{stat}_{dim}"
                    results[label] = stats[cond][stat][:, dim_index]


    def _combine_context_pvalues(self, results, test):
        sequence_context = self._config.get_sequence_context()
        sequence_context_weights = self._config.get_sequence_context_weights()
        if sequence_context_weights == "harmonic":
            # Generate weights as a symmetrical harmonic series
            weights = harmomic_series(sequence_context)
        else:
            weights = [1]*(2*sequence_context + 1)

        min_pos = results.pos.min()
        max_pos = results.pos.max()

        # We want to combine pvalues of neighbouring
        # positions, but the results dataframe may
        # have missing positions.
        # Additionally, applying a window function
        # means we lose the results for the positions
        # at the ends so we have to pad the results.
        padded_size = max_pos + 1 - min_pos + 2*sequence_context
        pvalues = np.ones(padded_size, dtype=float)
        indices = results.pos - min_pos
        pvalues[indices + sequence_context] = results[test]

        corr = cross_corr_matrix(pvalues[sequence_context:-sequence_context],
                                 sequence_context)

        window_size = 2*sequence_context + 1
        def func(window):
            if np.isnan(window[sequence_context]):
                return np.nan
            else:
                return combine_pvalues_hou(np.nan_to_num(window, nan=1),
                                           weights,
                                           corr)
        combinator = np.vectorize(func, signature='(n)->()')
        windows = sliding_window_view(pvalues, window_shape=window_size)
        combined_pvals = np.apply_along_axis(combinator, 1, windows)

        return combined_pvals[indices]


    def _get_cluster_counts(self, contingency, samples, conditions, predictions):
        depleted_cond = self._config.get_depleted_condition()
        depleted_cond_id = self._config.get_condition_ids()[depleted_cond]

        cond_counts = contingency.sum(2).to(torch.uint32)
        # contingency will be something like this:
        # e.g.      cluster
        # condition   0   1
        #         0  100  30
        #         1   90  40
        freq_contingency = contingency / cond_counts.unsqueeze(2)

        # We assume that the unmodified cluster is the one
        # in which the depleted condition has most of its
        # points. Hence, the modified is the one, where
        # the depleted condition has fewer points.
        # I.e. if the condition 0 is the depleted one,
        # in the example contingency above the modified
        # cluster will be cluster 1.
        mod_clusters = freq_contingency[:, depleted_cond_id, :].argmin(1)
        mod_clusters = mod_clusters.unsqueeze(1).expand_as(predictions)

        # Iterate all samples and calculate the
        # number of modified and non-modified reads.
        cluster_counts = {}
        for sample_label, sample in self._config.get_sample_ids().items():
            sample_predictions = predictions[:, samples == sample]
            sample_mod_clusters = mod_clusters[:, samples == sample]
            sample_totals = torch.sum(~sample_predictions.isnan(), dim=1)
            mod_counts = (sample_predictions == sample_mod_clusters).sum(dim=1)
            unmod_counts = sample_totals - mod_counts
            # unmod_count = (sample_predictions != sample_mod_clusters).sum().item()
            cluster_counts[f'{sample_label}_mod'] = mod_counts.cpu().numpy().astype(np.uint32)
            cluster_counts[f'{sample_label}_unmod'] = unmod_counts.cpu().numpy().astype(np.uint32)

        return cluster_counts


    def retry(self, fn, delay=5, backoff=5, max_attempts=3, exception=Exception):
        attempts = 0
        while attempts < max_attempts - 1:
            try:
                return fn()
            except exception as e:
                attempts += 1
                self._worker.log("warning", f"Got error {e} on attempt {attempts}/{max_attempts}")
                time.sleep(delay)
                delay += backoff
        else:
            return fn()


def nanstd(X, dim):
    nonnans = (~X.isnan()).sum(dim)
    return (((X - X.nanmean(dim).unsqueeze(1)) ** 2).nansum(dim) / nonnans) ** 0.5


def cross_corr_matrix(pvalues, context=2):
    """
    Calculate the cross correlation matrix of the
    pvalues for a given context.
    """
    if len(pvalues) < (3*context) + 3:
        raise RuntimeError("Not enough p-values for a context of order %s"%context)

    pvalues = np.nan_to_num(pvalues, nan=1)
    if any(pvalues == 0) or any(np.isinf(pvalues)) or any(pvalues > 1):
        raise RuntimeError("At least one p-value is invalid")

    if all(pvalues == 1):
        return(np.ones((2*context + 1, 2*context + 1)))

    matrix = np.zeros((2*context + 1, 2*context + 1))
    s = pvalues.size
    for i in range(-context, context + 1):
        for j in range(-context, i + 1):
            x = i + context
            y = j + context
            matrix[x, y] = np.corrcoef(
                    np.roll(pvalues, i)[context:s-context],
                    np.roll(pvalues, j)[context:s-context])[0, 1]
    return matrix + np.tril(matrix, -1).T


def harmomic_series(sequence_context):
    return [1/(abs(i) + 1)
            for i in range(-sequence_context, sequence_context + 1)]


def combine_pvalues_hou(pvalues, weights, cor_mat):
    """ Hou's method for the approximation for the distribution of the weighted
        combination of non-independent or independent probabilities.
        If any pvalue is nan, returns nan.
        https://doi.org/10.1016/j.spl.2004.11.028
        pvalues: list of pvalues to be combined
        weights: the weights of the pvalues
        cor_mat: a matrix containing the correlation coefficients between pvalues
        Test: when weights are equal and cor=0, hou is the same as Fisher
        print(combine_pvalues([0.1,0.02,0.1,0.02,0.3], method='fisher')[1])
        print(hou([0.1,0.02,0.1,0.02,0.3], [1,1,1,1,1], np.zeros((5,5))))
    """
    if len(pvalues) != len(weights):
        raise ValueError("Can't combine pvalues if pvalues and weights are not the same length.")
    if cor_mat.shape[0] != cor_mat.shape[1] or cor_mat.shape[0] != len(pvalues):
        raise ValueError("The correlation matrix needs to be square, with each" + \
                         "dimension equal to the length of the pvalued vector.")
    if all(p == 1 for p in pvalues):
        return 1
    if any(p == 0 or np.isinf(p) or p > 1 for p in pvalues):
        raise ValueError("At least one p-value is invalid")

    # Covariance estimation as in Kost and McDermott (eq:8)
    # https://doi.org/10.1016/S0167-7152(02)00310-3
    cov = lambda r: (3.263*r) + (0.710*r**2) + (0.027*r**3)
    k = len(pvalues)
    cov_sum = np.float64(0)
    sw_sum = np.float64(0)
    w_sum = np.float64(0)
    tau = np.float64(0)
    for i in range(k):
        for j in range(i + 1, k):
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
    combined_p_value = chi2.sf(tau/c, f)
    # Return a very small number if pvalue = 0
    if combined_p_value == 0:
        combined_p_value = np.finfo(np.float64).tiny
    return combined_p_value


def calculate_lors(contingencies):
    odds1 = (contingencies[:, 0, 0]/contingencies[:, 0, 1])
    odds2 = (contingencies[:, 1, 0]/contingencies[:, 1, 1])
    return torch.round(torch.log(odds1/odds2), decimals=3)


def get_contigency_matrices(conditions, predictions):
    num_positions = predictions.shape[0]
    contingencies = torch.zeros((num_positions, 2, 2))
    expanded_conditions = conditions.unsqueeze(0).expand_as(predictions)
    combs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for comb in combs:
        cond_val, pred_val = comb
        comb_count = torch.logical_and(
                expanded_conditions == cond_val,
                predictions == pred_val
            ).sum(1)
        contingencies[:, cond_val, pred_val] = comb_count
    return contingencies

