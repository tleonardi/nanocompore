#!/usr/bin/env python3
from __future__ import print_function

import pod5
import pysam
import numpy as np
import math
import pandas as pd
import itertools

from loguru import logger
from pkg_resources import resource_filename
from collections import Counter

from remora import io, refine_signal_map
from kmer import KmerData

from common import Kit, VAR_ORDER, NanocomporeError


RNA002_LEVELS_FILE = "kmer_models/rna002_5mer_levels_v1.txt"
RNA004_LEVELS_FILE = "kmer_models/rna004_9mer_levels_v1.txt"


class Remora:
    def __init__(self,
                 experiment,
                 config,
                 ref_id='',
                 start=0,
                 end=1,
                 seq='',
                 strand='+'):

        ########## Private fields ##########
        self._experiment = experiment
        self._config = config
        self._seq = seq
        self._min_coverage = config.get_min_coverage()
        self._ref_id = ref_id
        self._pod5_bam_tuples = []
        self._condition_labels = []
        self._sample_labels = []

        logger.trace(f"Creating Remora object for {ref_id}")
        self._build_pod5_bam_tuple()
        logger.trace(f"All pod5 and bam files were properly opened for {ref_id}")

        #Remora requires a kmer model file to resquiggle the data (Signal to sequence alignment)
        #This is defined as a singal refiner object in the Remora API
        #Without the signal refiner, it will not resquiggle, but instead merely return the ionic current stream
        try:
            self._sig_map_refiner = self._check_signal_refiner(kit=config.get_kit())
        except Exception as e:
            raise NanocomporeError ("failed to create the signal map refiner. Check that the kmer model table is up-to-date")

        #Remora requires a specific region in the reference to focus the resquiggling algorithm
        #This is part of the Remora API for creating the reference region
        try:
            self._ref_reg = io.RefRegion(ctg=ref_id, strand=strand, start=start, end=end)
        except:
            raise NanocomporeError (f"failed to create the remora reference region for {ref_id}")

        #TODO test Nanocompore accuracy using seconds per kmer or samples per kmer
        #Remora returns the number of datapoints per kmer, not the number of seconds the kmer persisted in the
        #sensitive region of the nanopore. This function uses the pod5 api to convert the sampling rate of the
        #sequencing (hz) to time per sample (seconds)
        try:
            self._time_per_sample = self._get_time_per_sample()
        except:
            raise NanocomporeError (f"failed to check for sampling rate. Likely something wrong with the pod5 file")


    def _get_samples_metrics(self):
        try:
            samples_metrics, bam_reads = self._remora_resquiggle(
                    self._config.get_downsample_high_coverage(),
                    self._config.get_kit())

            bam_reads = [[r.qname for r in reads] for reads in bam_reads]

            return samples_metrics, bam_reads
        except:
            raise NanocomporeError (f"failed to resquiggle with Remora")


    def _build_pod5_bam_tuple(self):
        for sample, pod5_fn, bam in self._experiment.get_sample_pod5_bam_data():
            try:
                pod5_fh = pod5.Reader(pod5_fn)
            except:
                raise NanocomporeError (f"failed to open pod5 file {pod5}")

            try:
                bam_fh = pysam.AlignmentFile(bam)
            except:
                raise NanocomporeError (f"failed to open bam file {bam}")

            self._pod5_bam_tuples.append((pod5_fh, bam_fh))
            self._sample_labels.append(sample)


    def _check_signal_refiner(self, kit):
        level_table = self._kmer_model_selector(kit=kit)
        self._kmer_size = self._kmer_size_detector(level_table=level_table)

        sig_map_refiner = refine_signal_map.SigMapRefiner(
            kmer_model_filename=level_table,
            scale_iters=0,
            algo='dwell_penalty',
            do_rough_rescale=True,
            do_fix_guage=True,
        )
        logger.debug(f"sig_map_refiner properly opened")
        return sig_map_refiner


    def _get_time_per_sample(self):
        logger.trace("Attempting to calculate time per sample")
        first_pod5 = self._pod5_bam_tuples[0][0]
        read = next(first_pod5.reads())
        sample_rate = read.run_info.sample_rate
        logger.trace(f"Sampling rate is {sample_rate} Hz")
        time_per_sample = 1.0/sample_rate
        logger.trace(f"Time per sample is {time_per_sample} seconds")
        return time_per_sample


    def _remora_resquiggle(self, max_reads, kit):
        RNA = 'RNA' in kit.name.upper()

        logger.info(f"Starting to resquiggle data with Remora API for {self._ref_reg.ctg}")
        samples_metrics, all_bam_reads = io.get_ref_reg_samples_metrics(
            self._ref_reg,
            self._pod5_bam_tuples,
            metric="dwell_trimmean_trimsd",
            sig_map_refiner=self._sig_map_refiner,
            max_reads=max_reads,
            reverse_signal=RNA,
            signal_type='norm',
        )
        logger.info(f"Data for {self._ref_reg.ctg} resquiggled")
        return samples_metrics, all_bam_reads


    def _kmer_model_selector(self, kit):
        if kit == Kit.RNA002:
            level_table = resource_filename('nanocompore', RNA002_LEVELS_FILE)
        elif kit == Kit.RNA004:
            level_table = resource_filename('nanocompore', RNA004_LEVELS_FILE)
        else:
            raise NotImplementedError(f"Kit {kit} not implemented yet.")

        return level_table


    def _kmer_size_detector(self, level_table):
        with open(level_table, 'r') as infile:
            line = infile.readline().strip().split('\t')
            kmer_size = len(line[0])
        logger.debug(f'The kmer size is {kmer_size}')
        return kmer_size


    def kmer_data_generator(self):
        # Resquiggle the signal and get the summary metrics for
        # each position of the transcript. The result is list of dicts
        # with one dict per sample and each dict has the type:
        # {metric: <numpy appray with shape (reads, positions)>, ...}
        samples_metrics, bam_reads = self._get_samples_metrics()

        # Get [(reads, positions, vars), ...] with one 3d tensor per sample
        per_sample_tensors = [np.stack([d[v] for v in VAR_ORDER], axis=2, dtype=np.float32)
                              for d in samples_metrics]
        # Delete the variables to allow GC to swipe unused
        # data while the generator is still used.
        del samples_metrics

        # Get [sample1_reads_count, sample2_reads_count, ...]
        sample_reads_count = [t.shape[0] for t in per_sample_tensors]

        # Get the label for each sample repeated by the number of reads
        sample_labels = np.array(self._sample_labels).repeat(sample_reads_count)

        # Get (samples*reads, pos, vars)
        tensor = np.concatenate(per_sample_tensors, axis=0, dtype=np.float32)
        del per_sample_tensors

        reads = None
        if self._config.get_read_level_data():
            reads = np.concatenate(bam_reads)

        for pos in range(self._ref_reg.len):
            kmer_seq = self._make_kmer_seq(pos)
            pos_data = tensor[:, pos, :]

            # Remove reads with nan values for any of the variables at that position
            non_nan_rows = ~np.isnan(pos_data).any(axis=1)
            pos_data = pos_data[non_nan_rows, :]
            pos_sample_labels = sample_labels[non_nan_rows]

            pos_reads = None
            if self._config.get_read_level_data():
                pos_reads = reads[non_nan_rows]

            kmer_data = KmerData(pos, kmer_seq, pos_data, pos_sample_labels, pos_reads, self._experiment)
            condition_counts = Counter(kmer_data.condition_labels)
            if all(condition_counts.get(cond, 0) >= self._min_coverage
                   for cond in self._experiment.get_condtion_labels()):
                yield kmer_data
            else:
                logger.trace(f'Skipping position {pos} due to insuffient coverage in both conditions')


    def _make_kmer_seq(self, pos):
        kmer = self._seq[pos:pos + self._kmer_size]
        kmer_diff = self._kmer_size - len(kmer)
        kmer = kmer + 'N'*kmer_diff
        if len(kmer) > self._kmer_size:
            raise NanocomporeError (f'{kmer} is longer than the modeled kmer size {self._kmer_size}')
        else:
            return kmer


    @property
    def kmer_size(self):
        return self._kmer_size

