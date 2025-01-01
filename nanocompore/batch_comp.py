import warnings
import gc

from collections import defaultdict

import time
import numpy as np
import pandas as pd
import statsmodels.discrete.discrete_model as dm
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

from nanocompore.gmm_statistics import fit_best_gmm


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
        t = time.time()
        data, samples, conditions, positions, transcripts = retry(lambda: self._kmers_to_tensor(kmers),
                                                                  exception=torch.OutOfMemoryError)

        # Standardize the data
        data = (data - data.nanmean(1).unsqueeze(1)) / self._nanstd(data, 1).unsqueeze(1)
        logger.debug(f"Finished converting kmers to tensor format ({time.time() - t})")

        auto_test_mask = None
        has_auto = 'auto' in self._config.get_comparison_methods()
        if has_auto:
            auto_test_mask = self._auto_test_mask(conditions)
        test_masks = self._get_test_masks(auto_test_mask, n_positions)

        results = pd.DataFrame({'transcript_id': transcripts,
                                'pos': positions})
        for test, mask in test_masks.items():
            logger.debug(f"Start {test}")
            t = time.time()
            test_results = self._run_test(test,
                                          data[mask, :, :],
                                          samples[mask, :],
                                          conditions[mask, :])
            logger.debug(f"Finished {test} ({time.time() - t})")
            self._merge_results(results, test_results, test, mask, auto_test_mask, n_positions)
        logger.debug(f"Start shift stats")
        t = time.time()
        retry(lambda: self._add_shift_stats(results, data, conditions),
              exception=torch.OutOfMemoryError)
        logger.debug(f"Finished shift stats ({time.time() - t})")

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

        device = self._config.get_device()
        data = torch.full((n_positions, n_reads, dims), np.nan, device=device)
        samples = torch.empty(n_positions, n_reads, dtype=torch.float32, device=device)
        conditions = torch.empty(n_positions, n_reads, dtype=torch.float32, device=device)
        positions = torch.empty(n_positions, dtype=int, device=device)
        transcripts = torch.empty(n_positions, dtype=int, device=device)

        for i, kmer in enumerate(kmers):
            motor_kmer = self._get_motor_kmer(kmer, kmers)
            reads = transcript_reads[kmer.transcript_id]
            read_order = np.vectorize(reads.get)(kmer.reads)
            data[i, :, INTENSITY] = self._order_values(n_reads, read_order, kmer.intensity)
            data[i, :, DWELL] = self._order_values(n_reads, read_order, kmer.dwell)
            if motor_kmer:
                motor_read_order = np.vectorize(reads.get)(motor_kmer.reads)
                data[i, :, MOTOR] = self._order_values(n_reads,
                                                       motor_read_order,
                                                       motor_kmer.dwell)
            samples[i, :] = self._order_values(n_reads, read_order, kmer.sample_ids)
            conditions[i, :] = self._order_values(n_reads, read_order, kmer.condition_ids)
            positions[i] = kmer.pos
            transcripts[i] = kmer.transcript_id
        return data.cpu(), samples.cpu(), conditions.cpu(), positions.cpu(), transcripts.cpu()


    def _order_values(self, n_reads, read_order, prop, dtype=torch.float32):
        device = self._config.get_device()
        values = torch.full((n_reads,), np.nan, device=device)
        values[read_order] = torch.tensor(prop, dtype=dtype, device=device)
        return values


    def _get_motor_kmer(self, kmer, kmers):
        for other_kmer in kmers:
            if (kmer.transcript_id == other_kmer.transcript_id and
                other_kmer.pos == kmer.pos + MOTOR_DWELL_EFFECT_OFFSET):
                return other_kmer
        return None


    def _get_test_masks(self, auto_test_mask, n_positions):
        base_tests = set(self._config.get_comparison_methods()) - {'auto'}
        additional_tests = set(auto_test_mask or []) - base_tests
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


    def _merge_results(self, results, test_results, test, mask, auto_test_mask, n_positions):
        has_auto = 'auto' in self._config.get_comparison_methods()

        if has_auto and 'auto_pvalue' not in results:
            results['auto_pvalue'] = np.full(n_positions, np.nan)
            results['auto_test'] = auto_test_mask
        
        for column, values in test_results.items():
            assigned_values = np.full(n_positions, np.nan)
            current_value = 0
            for i, included in enumerate(mask):
                if included:
                    assigned_values[i] = values[current_value]
                    current_value += 1
            results[column] = assigned_values
            if has_auto and 'pvalue' in column:
                results['auto_pvalue'] = np.where(auto_test_mask == test,
                                                  assigned_values,
                                                  results['auto_pvalue'])


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
            pval = gof_test_singlerep(test_data[i, :, :], conditions[i, :], self._config)
            pvals.append(pval)
        return {'GOF_pvalue': pvals}


    def _gof_test_multirep(self, test_data, samples, conditions):
        pvals = []
        for i in range(test_data.shape[0]):
            pval = gof_test_multirep(test_data[i, :, :], conditions[i, :], self._config)
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
        if dim3_data.shape[0] > 0:
            dim3_results = self._gmm_test_split(dim3_data,
                                                conditions[split == 0, :],
                                                indices[split == 0])
        if dim2_data.shape[0] > 0:
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
        if device == 'cuda':
            torch.cuda.empty_cache()
            s = test_data.element_size() * test_data.nelement()
            logger.info(f'Tensor size: {s}')
            logger.info(f'Memory reserved: {torch.cuda.memory_reserved()}')
            logger.info(f'Memory allocated: {torch.cuda.memory_allocated()}')
            logger.info(f'meminfo: {torch.cuda.mem_get_info(torch.device("cuda"))}')
            logger.info(f'X type: {test_data.dtype}')
            logger.info(f'X shape: {test_data.shape}')
        # s = time.time()
        # attempts = 0
        # exception = None
        # while attempts < 10:
        #     try:
        #         # gmm1 = GMM(n_components=1, device=device, reg_covar=1e-4)
        #         # gmm1.fit(test_data)
        #         # bic1 = gmm1.bic(test_data)
        #         # del gmm1
        #         # gc.collect()
        #         # if self._config.get_device() == 'cuda':
        #         #     torch.cuda.empty_cache()
        #         #     torch.cuda.reset_peak_memory_stats()
        #         gmm2 = GMM(n_components=2, device=device, reg_covar=1e-4)
        #         gmm2.fit(test_data)
        #         # bic2 = gmm2.bic(test_data)
        #         # gmm = GMM(n_components=2, device=device, reg_covar=1e-4, dtype=torch.bfloat16)
        #         # gmm.fit(test_data)
        #         break
        #     except torch.OutOfMemoryError as e:
        #         exception = e
        #         logger.info(f"GMM out of memory: {attempts}")
        #         attempts += 1
        #         time.sleep(5)
        # else:
        #     logger.info(f"GMM out of memory too many times: {attempts}")
        #     raise exception
        # logger.info(f"GMM fitting time: {time.time() - s}")
        # pred = gmm2.predict(test_data.to(torch.bfloat16))
        # del gmm2
        # gc.collect()
        # if self._config.get_device() == 'cuda':
        #     torch.cuda.empty_cache()
        #     torch.cuda.reset_peak_memory_stats()
        pvals = []
        lors = []
        dfs = 0
        gmms = 0
        logits = 0
        for i in range(test_data.shape[0]):

            s3 = time.time()

            valid = ~test_data[i, :, :].isnan().any(1)
            X = test_data[i, valid, :]
            try:
                gmm_fit = fit_best_gmm(X, max_components=2, cv_types=['full'], random_state=26)
                gmm_mod, gmm_type, gmm_ncomponents = gmm_fit
                gmms += time.time() - s3
            except Exception as e:
                logger.error(f"Fail in GMM fitting {e}")
                pvals.append(np.nan)
                lors.append(np.nan)
                continue
            
            if gmm_ncomponents < 2:
                pvals.append(1)
                lors.append(0)
                continue

            pred = gmm_mod.predict(X)

            # If the BIC with one component fits better
            # we conclude that there's no modification.
            # if bic1[i] <= bic2[i]:
            #     pvals.append(1)
            #     lors.append(0)
            #     continue

            # valid = ~conditions[i, :].isnan()
            s2 = time.time()

            # contingency = crosstab(conditions[i, valid], pred[i, valid]) + 1
            # dfs += time.time() - s2
            # pval = chi2_contingency(contingency).pvalue
            # pvals.append(pval)
            # lors.append(np.nan)

            X = np.ones((valid.sum() + 4, 2))
            X[:, 0] = np.append(conditions[i, valid], [0, 0, 1, 1])
            # y_pred = np.append(pred[i, valid].numpy(), [0, 0, 1, 1])
            y_pred = np.append(pred, [0, 0, 1, 1])
            dfs += time.time() - s2
            s = time.time()
            logit = dm.Logit(y_pred, X)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=PerfectSeparationWarning)
                warnings.filterwarnings('ignore', category=ConvergenceWarning)
                warnings.filterwarnings('ignore', category=RuntimeWarning)
                try:
                    # m = pd.DataFrame({'m': [','.join(xs) for xs in zip(pred[i, :], conditions[i, :])]}).m.value_counts()
                    # logger.info(f"matching: {m}")
                    logit_res = logit.fit(method='lbfgs', disp=False)
                    pval = logit_res.pvalues[0]
                    LOR = logit_res.params[0]
                    # pval = chi2_contingency(contingency_table).pvalue
                    # LOR = np.nan
                    pvals.append(pval)
                    lors.append(LOR)
                except Exception as e:
                    logger.error(f"Error in Logit: {e}")
                    logger.error(f"Error in Logit data: {y_pred} {y_pred.shape} {sum(y_pred)}")
                    pvals.append(np.nan)
                    lors.append(np.nan)
            logits += time.time() - s
        logger.info(f"GMM time: {gmms}")
        logger.info(f"Logit dataframe creation time: {dfs}")
        logger.info(f"Logit time: {logits}")
        return {'GMM_logit_pvalue': pvals,
                'GMM_logit_LOR': lors}


    def _split_by_ndim(self, test_data):
        any_motor_dwell_nan = (test_data[:, :, MOTOR].isnan() &
                               ~test_data[:, :, INTENSITY].isnan()).sum(1) > 0
        return (test_data[~any_motor_dwell_nan],
                test_data[any_motor_dwell_nan],
                any_motor_dwell_nan.int())


    def _add_shift_stats(self, results, data, conditions):
        data = data.to(self._config.get_device())

        stats = {0: {}, 1: {}}
        cond0_mask = (conditions == 0).to(self._config.get_device())
        cond0_data = data * cond0_mask.to(int).unsqueeze(2)
        stats[0]['mean'] = cond0_data.nanmean(1).to('cpu')
        stats[0]['median'] = cond0_data.nanmedian(1).values.to('cpu')
        stats[0]['std'] = self._nanstd(torch.where(cond0_mask.unsqueeze(2).repeat(1, 1, 3) == 0, data, np.nan), 1).to('cpu')

        cond1_mask = ~cond0_mask
        cond1_data = data * cond1_mask.to(int).unsqueeze(2)
        stats[1]['mean'] = cond1_data.nanmean(1).to('cpu')
        stats[1]['median'] = cond1_data.nanmedian(1).values.to('cpu')
        stats[1]['std'] = self._nanstd(torch.where(cond1_mask.unsqueeze(2).repeat(1, 1, 3) == 1, data, np.nan), 1).to('cpu')

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
        nonnans = ((~X.isnan()).to(int).sum(dim))
        return (((X - X.nanmean(dim).unsqueeze(1)) ** 2).nansum(dim) / nonnans) ** 0.5


def retry(fn, delay=5, backoff=5, max_attempts=10, exception=Exception):
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
        fn()


def crosstab_batch(a, b):
    B, N = a.shape
    table = torch.zeros((B, 2, 2), dtype=torch.int64)
    for i in range(2):
        for j in range(2):
            table[:, i, j] = ((a == i) & (b == j)).sum(dim=1)
    return table


def crosstab(a, b):
    table = torch.zeros((2, 2), dtype=torch.int64)
    for i in range(2):
        for j in range(2):
            table[i, j] = ((a == i) & (b == j)).sum()
    return table
    
