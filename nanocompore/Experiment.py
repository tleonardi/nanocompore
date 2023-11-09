# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import defaultdict
import logging
from loguru import logger
import random

# Third party
import numpy as np
import pandas as pd

# Local package
from nanocompore.common import *

class Experiment():
    def __init__(self, samples_tsv):
        self._input_data_df = self._build_input_data_df_from_tsv(samples_tsv)
        try:
            assert self._check_num_condition_labels()
        except:
            raise NanocomporeError(f"Incorrect number of condition labels. Nanocompore currently only supports 2 condition labels")

        try:
            assert self._check_unique_sample_labels
        except:
            raise NanocomporeError(f"One or more sample labels are not unique")

        self._build_condition_to_samples()
        self._build_sample_to_condition()

################### Private methods ###################

    #TODO fail fast with exception
    def _build_input_data_df_from_tsv(self, infile):
        if self._header_exists(infile):
            df = pd.read_csv(infile, sep='\t', header=0)
        else:
            df = pd.read_csv(infile, sep='\t', names=['Sample', 'Condition', 'pod5', 'bam'])
        return df

    def _header_exists(self, infile):
        with open(infile, 'r') as tsv:
            line = tsv.readline().strip().lower()
            if line == 'Sample\tCondition\tpod5\tbam':
                return True
            else:
                return False

    def _check_num_condition_labels(self):
        if len(set(self._input_data_df['Condition'].to_list())) == 2:
            return True
        else:
            return False
            
    def _check_unique_sample_labels(self):
        sample_labels = self._input_data_df['Sample'].to_list()
        if len(sample_labels) == len(set(sample_labels)):
            return True
        else:
            return False

    def _build_condition_to_samples(self):
        self._condition_to_samples = defaultdict(list)
        condition_labels = self._input_data_df['Condition'].to_list()
        sample_labels = self._input_data_df['Sample'].to_list()
        for condition, sample in zip(condition_labels, sample_labels):
            self._condition_to_samples[condition].append(sample)

    def _build_sample_to_condition(self):
        self._sample_to_condition = defaultdict()
        condition_labels = self._input_data_df['Condition'].to_list()
        sample_labels = self._input_data_df['Sample'].to_list()
        for condition, sample in zip(condition_labels, sample_labels):
            self._sample_to_condition[sample] = condition           

################### Public methods ###################

    def sample_to_condition(self, sample_label):
        return self._sample_to_condition[sample_label]

    def condition_to_samples(self, condition_label):
        return self._condition_to_samples[condition_label]

    def get_condtion_labels(self):
        return self._condition_to_samples.keys()

    def get_sample_labels(self):
        return self._sample_to_condition.keys()

    def check_bams_are_available(self):
        # Test that bam and bam.bai files are accessable
        logger.debug("Checking bam files can be accessed")
        for bam_fn in self._input_data_df['bam'].to_list():
            if not access_file(bam_fn):
                raise NanocomporeError("Cannot access bam file {}".format(bam_fn))

            bam_idx_fn = f"{bam_fn}.bai"
            if not access_file(bam_idx_fn):
                raise NanocomporeError("Cannot access bam index file {}".format(bam_fn))
        logger.info("All bams are accessable")
        return True

    def get_sample_pod5_bam_data(self):
        for sample, pod5, bam in self._get_input_data(data_labels=['Sample', 'pod5', 'bam']):
            yield sample, pod5, bam

    def get_sample_condition_bam_data(self):
        for sample, condition, bam in self._get_input_data(data_labels=['Sample', 'Condition', 'bam']):
            yield sample, condition, bam
        
    def _get_input_data(self, data_labels=['Sample', 'Condition', 'pod5', 'bam']):
        """
        Generator function that yields a tuple values in the order of data_labels input
        """
        for label in data_labels:
            if label not in ['Sample', 'Condition', 'pod5', 'bam']:
                raise NanocomporeError(f"{label} is not a support data label")

        for sample, condition, pod5, bam in self._input_data_df.itertuples(index=False):
            data = []
            for label in data_labels:
                if label == 'Sample':
                    data.append(sample)
                elif label == 'Condition':
                    data.append(condition)
                elif label == 'pod5':
                    data.append(pod5)
                elif label == 'bam':
                    data.append(bam)
            yield tuple(data)