# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import logging
from weakref import ref
from loguru import logger
import random

# Third party
import numpy as np
from tqdm import tqdm
from pyfaidx import Fasta
import pysam

# Local package
from nanocompore.common import *

# Set global random seed
downsample_random_seed = 42
random.seed(downsample_random_seed)

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist(object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
                 experiment,
                 config,
                 select_ref_id = [],
                 exclude_ref_id = []):
        """
        #########################################################
        * experiment
            object which contains methods relating to the sample and condition labels, 
            and the paths to the bam and pod5 files
        * config
            the input configuration object provided by the user
        * select_ref_id
            if given, only reference ids in the list will be selected for the analysis
        * exclude_ref_id
            if given, refid in the list will be excluded from the analysis
        """

        self._experiment = experiment
        self._fasta_fn = config.get_fasta_ref()


        # Filtering at ref_id level
        logger.info("Filtering out references with low coverage")
        self.ref_ids = self._select_ref(min_coverage=config.get_min_coverage(),
                                          min_ref_length=config.get_min_ref_length(),
                                          select_ref_ids=select_ref_id,
                                          exclude_ref_ids=exclude_ref_id)


    def __repr__(self):
        return "Whitelist: Number of references: {}".format(len(self))

    def __len__(self):
        return len(self.ref_ids)

    #~~~~~~~~~~~~~~PUBLIC METHODS AND PROPERTIES~~~~~~~~~~~~~~#
    @property
    def ref_id_list(self):
        return list(self.ref_ids)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _select_ref(self,
        min_coverage=10,
        min_ref_length=100,
        select_ref_ids=[],
        exclude_ref_ids=[]):
        """Select ref_id with a minimal coverage in both samples"""

        c = Counter()

        Condition_Counts, Sample_Counts = self._count_reads_in_bams()
        logger.info("Finished reading all the bams")

        valid_ref_ids = set()
        with Fasta(self._fasta_fn) as fasta:
            ref_ids = self._cross_check_ref_ids_with_select_and_exclude_ids(fasta.keys(),
                                                                            select_ref_ids,
                                                                            exclude_ref_ids)
            for ref_id in ref_ids:
                try:
                    # Discard reference transcripts shorter than the threshold
                    if len(fasta[ref_id]) < min_ref_length:
                        logger.trace(f"ref_id {ref_id} is shorter than min reference length: discarding it")
                        c["invalid_ref_id"] += 1
                        continue

                    # Discard reference transcripts that don't have enough coverage
                    conditions = Condition_Counts[ref_id].items()
                    if len(conditions) != 2:
                        logger.trace(f"ref_id {ref_id} doesn't have any coverage in both conditions: discarding it")
                        c["invalid_ref_id"] += 1
                        continue

                    if not all(count >= min_coverage for cond, count in conditions):
                        logger.trace(f"ref_id {ref_id} doesn't have any coverage in both conditions: discarding it")
                        c["invalid_ref_id"] += 1
                        continue

                    # If all valid add to set
                    logger.trace(f"ref_id {ref_id} has enough coverage in all samples: keeping it")
                    valid_ref_ids.add(ref_id)
                    # Save extra info for debug
                    c["valid_ref_id"] += 1

                    for sample in Sample_Counts[ref_id]:
                        c[sample] += Sample_Counts[ref_id][sample]

                except:
                    raise NanocomporeError(f"Error with testing coverage for {ref_id}")
                

        logger.debug(counter_to_str(c))
        logger.info("\tReferences remaining after reference coverage filtering: {}".format(len(valid_ref_ids)))
        return valid_ref_ids

    def _count_reads_in_bams(self):
        Condition_Counts = defaultdict(Counter)
        Sample_Counts = defaultdict(Counter)
        for sample, condition, bam_fn in self._experiment.get_sample_condition_bam_data():
            logger.info(f"Counting reference ids from {condition} {sample} bam file")
            with pysam.AlignmentFile(bam_fn, "rb") as bam:
                for read in bam:
                    Condition_Counts[read.reference_name][condition] += 1
                    Sample_Counts[read.reference_name][sample] += 1

        return Condition_Counts, Sample_Counts

    def _cross_check_ref_ids_with_select_and_exclude_ids(self, ref_ids, select_ref_ids=[], exclude_ref_ids=[]):
        valid_ids = set(ref_ids)

        if select_ref_ids or exclude_ref_ids:
            select_ref_ids = set(select_ref_ids)
            exclude_ref_ids = set(exclude_ref_ids)
            if select_ref_ids.intersection(exclude_ref_ids) == select_ref_ids.union(exclude_ref_ids):
                raise NanocomporeError("All of the reference ids to be selected are the same as the reference ids to be exluded")
            else:
                if select_ref_ids:
                    valid_ids = valid_ids.intersection(set(select_ref_ids))

                if exclude_ref_ids:
                    valid_ids = valid_ids.difference(set(exclude_ref_ids))
                
                if len(valid_ids) > 0:
                    return list(valid_ids)
                else:
                    raise NanocomporeError("There are no valid transcripts after cross checking with include and exlude parameters")
        else:
            return valid_ids