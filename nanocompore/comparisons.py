import warnings
import gc

from collections import defaultdict

import time
import numpy as np
import pandas as pd
import torch

from gmm_gpu.gmm import GMM
from loguru import logger
from scipy.stats import mannwhitneyu, ttest_ind, chi2_contingency
from scipy.stats.mstats import ks_twosamp
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

from nanocompore.common import MOTOR_DWELL_EFFECT_OFFSET
from nanocompore.common import NanocomporeError
from nanocompore.gof_tests import gof_test_multirep
from nanocompore.gof_tests import gof_test_singlerep


INTENSITY = 0
DWELL = 1
MOTOR = 2


class TranscriptComparator:
    """
    Compare a set of positions.
    """
    def __init__(self, config, random_seed=42, dtype=torch.float32):
        self._config = config
        self._random_seed = random_seed
        self._dtype = dtype


    def compare_transcript(self, transcript, kmers, device=None):
        if len(kmers) == 0:
            return (transcript, None)

        data, samples, conditions, positions = retry(lambda: self._kmers_to_tensor(kmers, device),
                                                     exception=torch.OutOfMemoryError)
        if device.startswith('cuda'):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()

        n_positions = len(positions)
        if n_positions == 0:
            return (transcript, None)

        results = pd.DataFrame({'transcript_id': transcript.id,
                                'pos': positions})

        logger.debug(f"Start shift stats")
        t = time.time()
        retry(lambda: self._add_shift_stats(results, data, conditions, device),
              exception=torch.OutOfMemoryError)
        if device.startswith('cuda'):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
        logger.debug(f"Finished shift stats ({time.time() - t})")

        # Standardize the data
        data = (data - data.nanmean(1).unsqueeze(1)) / self._nanstd(data, 1).unsqueeze(1)

        auto_test_mask = None
        has_auto = 'auto' in self._config.get_comparison_methods()
        if has_auto:
            auto_test_mask = self._auto_test_mask(conditions)
        test_masks = self._get_test_masks(auto_test_mask, n_positions)

        for test, mask in test_masks.items():
            logger.debug(f"Start {test}")
            t = time.time()
            test_results = self._run_test(test,
                                          data[mask, :, :],
                                          samples[mask, :],
                                          conditions[mask, :],
                                          device=device)
            logger.debug(f"Finished {test} ({time.time() - t})")
            self._merge_results(results, test_results, test, mask, auto_test_mask, n_positions)

        return transcript, results


    def _kmers_to_tensor(self, kmers, device):
        """
        Converts the list of kmers to tensors.
        Returns:
        - measurements: tensor with shape (positions, reads, vars) containing the
                        signal measurements.
        - samples: tensor with shape (positions, reads) containing sample ids.
        - conditions: tensor with shape (positions, reads) containing condition ids.
        - positions: tensor with shape (positions) with the transcript positions.
        """
        reads = {read
                 for kmer in kmers
                 for read in kmer.reads}

        min_pos, max_pos = np.inf, -1
        for kmer in kmers:
            if kmer.pos > max_pos:
                max_pos = kmer.pos
            if kmer.pos < min_pos:
                min_pos = kmer.pos

        read_indices = dict(zip(reads, range(len(reads))))
        get_indices = np.vectorize(read_indices.get)

        initial_positions = max_pos - min_pos + 1

        tensor = np.full((initial_positions, len(reads), 6), np.nan)
        valid_positions = torch.full((initial_positions,), False)
        tmp_matrix = np.empty((len(reads), 5), dtype=float)
        for kmer in sorted(kmers, key=lambda kmer: kmer.pos, reverse=True):
            indices = get_indices(kmer.reads)
            pos = kmer.pos - min_pos
            valid_positions[pos] = True
            nreads = len(kmer.reads)
            # It's a bit faster to put all the kmer data
            # (which is numpy arrays), into one big matrix
            # in consecutive positions and then to reorder
            # the rows appropritately.
            tmp_matrix[:nreads, 0] = kmer.reads
            tmp_matrix[:nreads, 1] = kmer.condition_ids
            tmp_matrix[:nreads, 2] = kmer.sample_ids
            tmp_matrix[:nreads, 3] = kmer.intensity
            tmp_matrix[:nreads, 4] = kmer.dwell
            tensor[pos, indices, :5] = tmp_matrix[:nreads, :]
        end = initial_positions - MOTOR_DWELL_EFFECT_OFFSET
        if end > 0:
            tensor[:end, :, 5] = tensor[MOTOR_DWELL_EFFECT_OFFSET:initial_positions, :, 4]

        # valid positions are those that have at least some values
        tensor = torch.tensor(tensor[valid_positions.numpy()],
                              dtype=self._dtype)

        positions = torch.arange(min_pos, max_pos + 1)[valid_positions]
        return (tensor[:, :, 3:6].clone(), # measurements
                tensor[:, :, 2].clone(), # samples
                tensor[:, :, 1].clone(), # conditions
                positions.clone())


    def _get_motor_kmer(self, kmer, kmers):
        for other_kmer in kmers:
            if (kmer.transcript_id == other_kmer.transcript_id and
                other_kmer.pos == kmer.pos + MOTOR_DWELL_EFFECT_OFFSET):
                return other_kmer
        return None


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


    def _auto_test_mask(self, conditions):
        depleted_reads = (conditions == 0).nansum(1)
        non_depleted_reads = (conditions == 1).nansum(1)
        return np.array([self._resolve_auto_test(counts)
                         for counts in zip(depleted_reads, non_depleted_reads)])


    def _resolve_auto_test(self, condition_counts):
        min_coverage = min(condition_counts)
        if min_coverage < 256:
            return 'KS'
        else:
            return 'GMM'


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
            return self._gmm_test(test_data, conditions, device)
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
            stat_test = lambda x,y: mannwhitneyu(x, y, alternative='two-sided')
        elif test in ["kolmogorov_smirnov", "KS"]:
            stat_test = ks_twosamp
        elif test in ["t_test", "TT"]:
            stat_test = lambda x,y: ttest_ind(x, y, equal_var=False)
        else:
            raise NanocomporeError("Invalid statistical method name (MW, KS, TT)")

        cond0_mask = conditions == 0

        intensity_pvals = []
        dwell_pvals = []
        for i in range(test_data.shape[0]):
            cond0_data = test_data[i, cond0_mask[i], :]
            cond1_data = test_data[i, ~cond0_mask[i], :]
            try:
                pval = stat_test(self._drop_nans(cond0_data[:, INTENSITY]),
                                 self._drop_nans(cond1_data[:, INTENSITY])).pvalue
                intensity_pvals.append(pval)
            except:
                intensity_pvals.append(np.nan)
            try:
                pval = stat_test(self._drop_nans(cond0_data[:, DWELL]),
                                 self._drop_nans(cond1_data[:, DWELL])).pvalue
                dwell_pvals.append(pval)
            except:
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
            pval = gof_test_singlerep(test_data[i, :, :], conditions[i, :], self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gof_test_multirep(self, test_data, samples, conditions):
        pvals = []
        for i in range(test_data.shape[0]):
            pval = gof_test_multirep(test_data[i, :, :], conditions[i, :], self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gmm_test(self, test_data, conditions, device):
        # For each position some reads have motor dwell
        # and others don't. We split the positions into
        # two sets: one set of positions where all reads
        # have motor dwell time - those would be tested
        # with a 3D GMM - and a second set of positions
        # where some reads are missing the motor dwell -
        # those would be tested with a 2D GMM.
        dim3_data, dim2_data, split = self._split_by_ndim(test_data)
        indices = torch.arange(test_data.shape[0])
        pvals = []
        lors = []
        if dim3_data.shape[0] > 0:
            dim3_results = self._gmm_test_split(dim3_data,
                                                conditions[split == 0, :],
                                                indices[split == 0],
                                                device)
        if dim2_data.shape[0] > 0:
            dim2_results = self._gmm_test_split(dim2_data,
                                                conditions[split == 1, :],
                                                indices[split == 1],
                                                device)
        dim3_i = 0
        dim2_i = 0
        for s in split:
            if s == 0:
                pvals.append(dim3_results['GMM_chi2_pvalue'][dim3_i])
                lors.append(dim3_results['GMM_LOR'][dim3_i])
                dim3_i += 1
            else:
                pvals.append(dim2_results['GMM_chi2_pvalue'][dim2_i])
                lors.append(dim2_results['GMM_LOR'][dim2_i])
                dim2_i += 1
        return {'GMM_chi2_pvalue': pvals,
                'GMM_LOR': lors}


    def _gmm_test_split(self, data, conditions, indices, device):
        test_data = data.to(self._dtype)
        s = time.time()
        def fit_model(components):
            gmm = GMM(n_components=components,
                      device=device,
                      random_seed=self._random_seed,
                      dtype=self._dtype)
            gmm.fit(test_data)
            return gmm
        gmm1 = retry(lambda: fit_model(1), exception=torch.OutOfMemoryError)
        bic1 = gmm1.bic(test_data)
        del gmm1
        gc.collect()
        gmm2 = retry(lambda: fit_model(2), exception=torch.OutOfMemoryError)

        logger.info(f"GMM fitting time: {time.time() - s}")
        pred = gmm2.predict(test_data)
        bic2 = gmm2.bic(test_data)
        del gmm2
        gc.collect()
        if device.startswith('cuda'):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
        lors = self._calculate_log_odds_ratios(pred, conditions)
        pvals = []
        for i in range(test_data.shape[0]):
            if bic1[i] <= bic2[i]:
                pvals.append(np.nan)
                lors[i] = np.nan
                continue
            valid = ~conditions[i, :].isnan()
            contingency = crosstab(conditions[i, valid], pred[i, valid]) + 1
            pval = chi2_contingency(contingency).pvalue
            pvals.append(pval)
        return {'GMM_chi2_pvalue': pvals,
                'GMM_LOR': lors}


    def _calculate_log_odds_ratios(self, pred, conditions):
        cond0 = (conditions == 0).any(0)
        reads00_num = (pred[:, cond0] == 0).sum(1) + 1
        reads01_num = pred[:, cond0].nansum(1) + 1
        cond1 = (conditions == 1).any(0)
        reads10_num = (pred[:, cond1] == 0).sum(1) + 1
        reads11_num = pred[:, cond1].nansum(1) + 1
        return torch.log((reads00_num/reads01_num)/(reads10_num/reads11_num))


    def _split_by_ndim(self, test_data):
        any_motor_dwell_nan = (test_data[:, :, MOTOR].isnan() &
                               ~test_data[:, :, INTENSITY].isnan()).sum(1) > 0
        return (test_data[~any_motor_dwell_nan],
                test_data[any_motor_dwell_nan, :, :MOTOR],
                any_motor_dwell_nan.int())


    def _add_shift_stats(self, results, data, conditions, device):
        data = data.to(device)

        stats = {0: {}, 1: {}}
        cond0_mask = (conditions == 0).unsqueeze(2).to(device)

        cond0_data = torch.where(cond0_mask, data, np.nan)
        stats[0]['mean'] = cond0_data.nanmean(1).cpu()
        stats[0]['median'] = cond0_data.nanmedian(1).values.cpu()
        stats[0]['std'] = self._nanstd(cond0_data, 1).cpu()

        cond1_data = torch.where(~cond0_mask, data, np.nan)
        stats[1]['mean'] = cond1_data.nanmean(1).cpu()
        stats[1]['median'] = cond1_data.nanmedian(1).values.cpu()
        stats[1]['std'] = self._nanstd(cond1_data, 1).cpu()

        dims = {
                'intensity': INTENSITY,
                'dwell': DWELL,
                'motor_dwell': MOTOR
               }

        for cond in [0, 1]:
            for stat in ['mean', 'median', 'std']:
                for dim, dim_index in dims.items():
                    label = f"c{cond+1}_{stat}_{dim}"
                    results[label] = stats[cond][stat][:, dim_index]


    def _nanstd(self, X, dim):
        nonnans = (~X.isnan()).to(int).sum(dim)
        return (((X - X.nanmean(dim).unsqueeze(1)) ** 2).nansum(dim) / nonnans) ** 0.5


def retry(fn, delay=5, backoff=5, max_attempts=3, exception=Exception):
    attempts = 0
    while attempts < max_attempts - 1:
        try:
            return fn()
        except exception as e:
            attempts += 1
            logger.warning(f"Got error {e} on attempt {attempts}/{max_attempts}")
            time.sleep(delay)
            delay += backoff
    else:
        return fn()


def crosstab(a, b):
    table = torch.zeros((2, 2), dtype=torch.int)
    for i in range(2):
        for j in range(2):
            table[i, j] = ((a == i) & (b == j)).sum()
    return table

