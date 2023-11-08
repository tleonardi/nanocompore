import itertools
from lib2to3.pgen2.token import RPAR
import nanocompore.Sample as Sample
from loguru import logger
from nanocompore.common import NanocomporeError
import numpy as np
import pandas as pd
from collections import Counter, defaultdict


import nanocompore.Remora as Remora


#TODO change to just "Transcript"
class Transcript():
    def __init__(self, ref_id='',
                 experiment = '',
                 ref_seq=0,
                 min_coverage=30,
                 max_coverage=float('NaN'),
                 kit='RNA002'):

        logger.trace(f"Creating transcript object for {ref_id}")
        self._experiment = experiment
        self._name = ref_id
        self._ref_seq = ref_seq
        self._length = len(self._ref_seq)

        self._remora_object = Remora.Remora(experiment=self._experiment, ref_id=ref_id, start=0, end=self._length, seq=ref_seq, 
                                            min_reads=min_coverage, max_reads=max_coverage, strand='+', kit=kit)

    ########## Public ##########
    @property
    def name(self):
        return self._name
    @property
    def length(self):
        return self._length
    @property
    def seq(self):
        return self._ref_seq
    @property
    def condition_labels(self):
        return list(self._experiment.get_condtion_labels())
    @property
    def sample_labels(self):
        return self._experiment.get_sample_labels()

    def sample_2_condition(self, sample_label):
        return self._experiment.sample_to_condition(sample_label)

    def generate_data(self):
        logger.trace(f"Iterating through kmer data positions for {self._name}")
        for kmer_data in self._remora_object.kmer_data_generator():
            yield kmer_data