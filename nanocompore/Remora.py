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

from remora import io, refine_signal_map

from nanocompore.common import *

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

        try:
            self._samples_metrics = self._remora_resquiggle(config.get_downsample_high_coverage(), config.get_kit())
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

    def _remora_resquiggle(self, max_reads, kit='RNA002'):
        if 'RNA' in kit.upper():
            RNA = True
        else:
            RNA = False
        
        logger.info(f"Starting to resquiggle data with Remora API for {self._ref_reg.ctg}")
        samples_metrics, all_bam_reads = io.get_ref_reg_samples_metrics(
            self._ref_reg,
            self._pod5_bam_tuples,
            metric="dwell_trimmean_trimsd",
            sig_map_refiner=self._sig_map_refiner,
            max_reads=max_reads,
            reverse_signal=RNA
        )
        logger.info(f"Data resquiggled")
        return samples_metrics

    def _kmer_model_selector(self, kit=('RNA002')):#, 'RNA004', 'DNA260', 'DNA400')):
        if kit not in ('RNA002'):#, 'RNA004', 'DNA260', 'DNA400')
            raise NanocomporeError (f'Unsupported kit, exiting')

        #TODO make the kmer_models path work
        else:
            if kit == 'RNA002':
                level_table = resource_filename("nanocompore", "rna002_5mer_levels_v1.txt")
                #level_table = resource_filename("nanocompore/kmer_models/rna_r9.4_180mv_70bps", "5mer_levels_v1.txt")
            elif kit == 'RNA004':
                level_table = resource_filename("nanocompore/kmer_models/rna004", "9mer_levels_v1.txt")
            elif kit == 'DNA260':
                level_table = resource_filename("nanocompore/kmer_models/dna_r10.4.1_e8.2_260bps", "9mer_levels_v1.txt")
            elif kit == 'DNA400':
                level_table = resource_filename("nanocompore/kmer_models/dna_r10.4.1_e8.2_400bps", "9mer_levels_v1.txt")

        return level_table

    def _kmer_size_detector(self, level_table):
        with open(level_table, 'r') as infile:
            line = infile.readline().strip().split('\t')
            kmer_size = len(line[0])
        logger.debug(f'The kmer size is {kmer_size}')
        return kmer_size

    def kmer_data_generator(self):
        logger.debug(f"Iterating through each position of {self._ref_reg.ctg} and return a dataframe for each kmer raw ionic current data")
        for pos in range(self._ref_reg.len):
            kmer_seq = self._make_kmer_seq(pos)
            intensity = []
            dwell = []
            sample_labels = []
            condition_labels = []
            for i, sample in enumerate(self._samples_metrics):
                sample_label = self._sample_labels[i]
                condition_label = self._experiment.sample_to_condition(sample_label)

                intensity.append(sample['trimmean'][:, pos])
                dwell.append(sample['dwell'][:, pos])

                sample_labels.append([sample_label]*len(sample['trimmean'][:, pos]))
                condition_labels.append([condition_label]*len(sample['trimmean'][:, pos]))

            intensity = np.array(self._unpack_list(intensity))
            dwell = np.array(self._unpack_list(dwell)) * self._time_per_sample
            sample_labels = self._unpack_list(sample_labels)
            condition_labels = self._unpack_list(condition_labels)

            dataframe = self._build_dataframe(intensity, dwell, sample_labels, condition_labels)
            if self._check_condition_counts_threshold(dataframe):
                kmer_data = Remora._Kmer_Data(pos, kmer_seq, dataframe)
                yield kmer_data
            else:
                logger.trace(f'Skipping position {pos} due to insuffient coverage in both conditions')

    def _build_dataframe(self, intensity, dwell, sample_labels, condition_labels):
        dataframe = pd.DataFrame([intensity, dwell, sample_labels, condition_labels]).transpose()
        dataframe.columns = ('intensity', 'dwell', 'sample', 'condition')
        dataframe.dropna(inplace=True)

        dataframe.dwell = dataframe.dwell.apply(math.log10)
        return dataframe

    def _unpack_list(self, list_of_lists):
        return list(itertools.chain(*list_of_lists))

    def _check_condition_counts_threshold(self, data_frame, condition='condition'):
        condition_counts = data_frame[condition].value_counts()
        return all(count >= self._min_coverage for count in condition_counts)  

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


    class _Kmer_Data():
        def __init__(self, pos, kmer, data):
            self._kmer = kmer
            self._data = data
            self._pos = pos

        ########## Public ##########
        @property
        def pos(self):
            return self._pos
        @property
        def kmer(self):
            return self._kmer
        @property
        def intensity(self):
            return self._data['intensity'].to_numpy()
        @property
        def dwell(self):
            return self._data['dwell'].to_numpy()
        @property
        def sample_labels(self):
            return self._data['sample'].tolist()
        @property
        def condition_labels(self):
            return self._data['condition'].tolist()

        def get_condition_kmer_data(self, condition_label, data_type=''):
            if data_type not in ('intensity', 'dwell'):
                raise NanocomporeError (f"{data_type} is not currently supported")

            return self._data[self._data['condition'] == condition_label][data_type].to_numpy()

        def get_condition_kmer_intensity_data(self, condition_label):
            return self.get_condition_kmer_data(condition_label, data_type='intensity')

        def get_condition_kmer_dwell_data(self, condition_label):
            return self.get_condition_kmer_data(condition_label, data_type='dwell')