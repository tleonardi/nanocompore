import warnings

from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.discrete.discrete_model as dm
import torch

from gmm_gpu.gmm import GMM
from loguru import logger
from scipy.stats import mannwhitneyu, ttest_ind
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


class BatchComp:
    """
    Compare a set of positions.
    """
    def __init__(self, config):
        self._config = config


    def compare_transcripts(self, kmers):
        n_positions = len(kmers)
        logger.debug("Start converting kmers to tensor format")
        data, samples, conditions, positions, transcripts = self._kmers_to_tensor(kmers)
        logger.debug("Finished converting kmers to tensor format")

        auto_test_mask = self._auto_test_mask(conditions)
        test_masks = self._get_test_masks(auto_test_mask, n_positions)

        results = pd.DataFrame({'transcript_id': transcripts,
                                'pos': positions})
        for test, mask in test_masks.items():
            test_results = self._run_test(test,
                                          data[mask, :, :],
                                          samples[mask, :],
                                          conditions[mask, :])
            self._merge_results(results, test_results, mask, n_positions)
        self._add_shift_stats(results, data, conditions)

        return results


    def _kmers_to_tensor(self, kmers):
        # transcript_reads has type:
        # transcript_id -> {read_id -> index}
        # We'll use it to ensure the read order is the
        # same for all positions of the same transcript.
        transcript_reads = self._transcript_reads(kmers)

        n_positions = len(kmers)
        if len(transcript_reads) > 0:
            n_reads = max(len(reads) for reads in transcript_reads.values())
        else:
            n_reads = 0
        dims = 3

        data = torch.full((n_positions, n_reads, dims), np.nan)
        samples = torch.empty(n_positions, n_reads, dtype=int)
        conditions = torch.empty(n_positions, n_reads, dtype=int)
        positions = torch.empty(n_positions, dtype=int)
        transcripts = torch.empty(n_positions, dtype=int)

        for i, kmer in enumerate(kmers):
            motor_kmer = self._get_motor_kmer(kmer, kmers)
            reads = transcript_reads[kmer.transcript_id]
            read_order = np.vectorize(reads.get)(kmer.reads)
            data[i, :, INTENSITY] = self._order_values(n_reads, read_order, kmer, kmer.intensity)
            data[i, :, DWELL] = self._order_values(n_reads, read_order, kmer, kmer.dwell)
            if motor_kmer:
                motor_read_order = np.vectorize(reads.get)(motor_kmer.reads)
                data[i, :, MOTOR] = self._order_values(n_reads,
                                                       motor_read_order,
                                                       motor_kmer,
                                                       motor_kmer.dwell)
            samples[i, :] = self._order_values(n_reads, read_order, kmer, kmer.sample_ids)
            conditions[i, :] = self._order_values(n_reads, read_order, kmer, kmer.condition_ids)
            positions[i] = kmer.pos
            transcripts[i] = kmer.transcript_id
        return data, samples, conditions, positions, transcripts


    def _order_values(self, n_reads, read_order, kmer, prop):
        values = torch.full((n_reads,), np.nan)
        values[read_order] = torch.tensor(prop, dtype=torch.float32)
        return values


    def _get_motor_kmer(self, kmer, kmers):
        for other_kmer in kmers:
            if (kmer.transcript_id == other_kmer.transcript_id and
                other_kmer.pos == kmer.pos + MOTOR_DWELL_EFFECT_OFFSET):
                return other_kmer
        return None


    def _get_test_masks(self, auto_test_mask, n_positions):
        base_tests = set(self._config.get_comparison_methods()) - {'auto'}
        additional_tests = set(auto_test_mask) - base_tests
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
        elif min_coverage < 1024:
            return 'GMM'
        elif min_coverage >= 1024:
            return 'GOF'
        return 'KS'


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


    def _run_test(self, test, test_data, samples, conditions):
        if test.upper() in ['MW', 'KS', 'TT']:
            return self._nonparametric_test(test, test_data, conditions)
        elif test.upper() == 'GMM':
            return self._gmm_test(test_data, conditions)
        elif test.upper() == 'GOF':
            return self._gof_test(test_data, samples, conditions)
        else:
            raise NotImplementedError(f"Test {test} is not implemented!")


    def _merge_results(self, results, test_results, mask, n_positions):
        for column, values in test_results.items():
            assigned_values = np.full(n_positions, np.nan)
            i = 0
            for included in mask:
                if included:
                    assigned_values[i] = values[i]
                    i += 1
            results[column] = assigned_values


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
            pval = stat_test(self._drop_nans(cond0_data[:, INTENSITY]),
                             self._drop_nans(cond1_data[:, INTENSITY])).pvalue
            intensity_pvals.append(pval)
            pval = stat_test(self._drop_nans(cond0_data[:, DWELL]),
                             self._drop_nans(cond1_data[:, DWELL])).pvalue
            dwell_pvals.append(pval)
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
            pval = gof_test_multirep(test_data[i, :, :], conditions, self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gmm_test(self, test_data, conditions):
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
        dim3_results = self._gmm_test_split(dim3_data,
                                            conditions[split == 0, :],
                                            indices[split == 0])
        dim2_results = self._gmm_test_split(dim2_data,
                                            conditions[split == 1, :],
                                            indices[split == 1])
        dim3_i = 0
        dim2_i = 0
        for s in split:
            if s == 0:
                pvals.append(dim3_results['GMM_logit_pvalue'][dim3_i])
                lors.append(dim3_results['GMM_logit_LOR'][dim3_i])
                dim3_i += 1
            else:
                pvals.append(dim2_results['GMM_logit_pvalue'][dim2_i])
                lors.append(dim2_results['GMM_logit_LOR'][dim2_i])
                dim2_i += 1
        return {'GMM_logit_pvalue': pvals,
                'GMM_logit_LOR': lors}


    def _gmm_test_split(self, test_data, conditions, indices):
        device = self._config.get_device()
        gmm = GMM(n_components=2, device=device, reg_covar=1e-4)
        gmm.fit(test_data)
        pred = gmm.predict(test_data)
        pvals = []
        lors = []
        for i in range(test_data.shape[0]):
            X = pd.DataFrame({'component': conditions[i, :],
                              'intercept': 1})
            logit = dm.Logit(pred[i, :], X)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=PerfectSeparationWarning)
                warnings.filterwarnings('ignore', category=ConvergenceWarning)
                try:
                    logit_res = logit.fit(disp=False)
                    pval = logit_res.pvalues.iloc[1]
                    LOR = logit_res.params.iloc[1]
                    pvals.append(pval)
                    lors.append(LOR)
                except ConvergenceWarning:
                    logger.warn('Caught ConvergenceWarning')
                    pvals.append(np.nan)
                    lors.append(np.nan)
        return {'GMM_logit_pvalue': pvals,
                'GMM_logit_LOR': lors}


    def _split_by_ndim(self, test_data):
        any_motor_dwell_nan = (test_data[:, :, MOTOR].isnan() &
                               ~test_data[:, :, INTENSITY].isnan()).sum(1) > 0
        return (test_data[~any_motor_dwell_nan],
                test_data[any_motor_dwell_nan],
                any_motor_dwell_nan.int())


    def _add_shift_stats(self, results, data, conditions):
        cond0_mask = conditions == 0
        stats = {
                'mean': torch.nanmean,
                'median': torch.nanmedian,
                'std': lambda v: self._drop_nans(v).std()
                }
        dims = {
                'intensity': INTENSITY,
                'dwell': DWELL,
                'motor_dwell': MOTOR
               }

        for i in range(data.shape[0]):
            cond0_data = data[:, cond0_mask[i], :]
            cond1_data = data[:, ~cond0_mask[i], :]

            conds = {
                    'c1': cond0_data,
                    'c2': cond1_data
                    }

            for cond in conds:
                for stat in stats:
                    for dim in dims:
                        value = stats[stat](conds[cond][:, dims[dim]])
                        label = f"{cond}_{stat}_{dim}"
                        if label not in results.columns:
                            results[label] = np.full(data.shape[0], np.nan)
                        results.loc[i, label] = value.item()

