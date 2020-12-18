# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import logging
from loguru import logger
import random

# Third party
import numpy as np
from tqdm import tqdm
from pyfaidx import Fasta

# Local package
from nanocompore.common import *

# Set global random seed
downsample_random_seed = 42
random.seed(downsample_random_seed)

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist(object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
        eventalign_fn_dict,
        fasta_fn,
        min_coverage = 10,
        min_ref_length = 100,
        downsample_high_coverage = False,
        max_invalid_kmers_freq = 0.1,
        max_NNNNN_freq = 0.1,
        max_mismatching_freq = 0.1,
        max_missing_freq = 0.1,
        select_ref_id = [],
        exclude_ref_id = []):
        """
        #########################################################
        * eventalign_fn_dict
            Multilevel dictionnary indicating the condition_label, sample_label and file name of the eventalign_collapse output
            example d = {"S1": {"R1":"path1.tsv", "R2":"path2.tsv"}, "S2": {"R1":"path3.tsv", "R2":"path4.tsv"}}
        * fasta_fn
            Path to a fasta file corresponding to the reference used for read alignemnt
        * min_coverage
            minimal coverage required in both samples
        * min_ref_length
            minimal length of a reference transcript to be considered in the analysis
        * downsample_high_coverage
            For reference with higher coverage, downsample by randomly selecting reads.
        * max_invalid_kmers_freq
            maximum frequency of NNNNN, mismatching and missing kmers in reads
            If None, then the max_NNNNN_freq, max_mismatching_freq, max_missing_freq agrs will be used instead
        * max_NNNNN_freq
            maximum frequency of NNNNN kmers in reads (1 to deactivate)
        * max_mismatching_freq
            maximum frequency of mismatching kmers in reads (1 to deactivate)
        * max_missing_freq
            maximum frequency of missing kmers in reads (1 to deactivate)
        * select_ref_id
            if given, only reference ids in the list will be selected for the analysis
        * exclude_ref_id
            if given, refid in the list will be excluded from the analysis
        """

        # Check index files
        self.__filter_invalid_kmers = True
        for sample_dict in eventalign_fn_dict.values():
            for fn in sample_dict.values():
                idx_fn = fn+".idx"
                if not access_file(idx_fn):
                    raise NanocomporeError("Cannot access eventalign_collapse index file {}".format(idx_fn))
                # Check header line and set a flag to skip filter if the index file does not contain kmer status information
                with open(idx_fn, "r") as fp:
                    header = fp.readline().rstrip().split("\t")
                if not all_values_in (("ref_id", "read_id", "byte_offset", "byte_len"), header):
                    raise NanocomporeError("The index file {} does not contain the require header fields".format(idx_fn))
                if not all_values_in (("kmers", "NNNNN_kmers", "mismatch_kmers", "missing_kmers"), header):
                    self.__filter_invalid_kmers = False
                    logger.debug("Invalid kmer information not available in index file")

        self.__eventalign_fn_dict = eventalign_fn_dict

        # Get number of samples
        n = 0
        for sample_dict in self.__eventalign_fn_dict.values():
            for sample_lab in sample_dict.keys():
                n+=1
        self.__n_samples = n

        # Test is Fasta can be opened
        try:
            with Fasta(fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        # Create reference index for both files
        logger.info("Reading eventalign index files")
        ref_reads = self.__read_eventalign_index(
            eventalign_fn_dict = eventalign_fn_dict,
            max_invalid_kmers_freq = max_invalid_kmers_freq,
            max_NNNNN_freq = max_NNNNN_freq,
            max_mismatching_freq = max_mismatching_freq,
            max_missing_freq = max_missing_freq,
            select_ref_id = select_ref_id,
            exclude_ref_id = exclude_ref_id)

        # Filtering at transcript level
        logger.info("Filtering out references with low coverage")
        self.ref_reads = self.__select_ref(
            ref_reads = ref_reads,
            min_coverage=min_coverage,
            min_ref_length=min_ref_length,
            downsample_high_coverage=downsample_high_coverage)

        self.__min_coverage = min_coverage
        self.__min_ref_length = min_ref_length
        self.__downsample_high_coverage = downsample_high_coverage
        self.__max_invalid_kmers_freq = max_invalid_kmers_freq

    def __repr__(self):
        return "Whitelist: Number of references: {}".format(len(self))

    def __str__(self):
        m = ""
        for ref_id, ref_dict in self.ref_reads.items():
            m += "{}\n".format(ref_id)
            for cond_lab, cond_dict in ref_dict.items():
                for sample_lab, read_list in cond_dict.items():
                    m += "\t{} {}\n".format(cond_lab, sample_lab)
                    for reads in read_list:
                        m += "\t\t"
                        for k,v in reads.items():
                            m += "{}:{}  ".format(k,v)
                        m += "\n"
        return m

    def __len__(self):
        return len(self.ref_reads)

    def __iter__(self):
        for i, j in self.ref_reads.items():
            yield(i,j)

    def __getitem__(self, items):
        return self.ref_reads.get(items, None)

    #~~~~~~~~~~~~~~PUBLIC METHODS AND PROPERTIES~~~~~~~~~~~~~~#
    @property
    def ref_id_list(self):
        return list(self.ref_reads.keys())

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def __read_eventalign_index(self,
        eventalign_fn_dict,
        max_invalid_kmers_freq,
        max_NNNNN_freq,
        max_mismatching_freq,
        max_missing_freq,
        select_ref_id,
        exclude_ref_id):
        """Read the 2 index files and sort by sample and ref_id in a multi level dict"""

        ref_reads = OrderedDict()

        for cond_lab, sample_dict in eventalign_fn_dict.items():
            for sample_lab, fn in sample_dict.items():
                idx_fn = fn+".idx"
                with open(idx_fn) as fp:

                    # Get column names from header
                    col_names = fp.readline().rstrip().split()
                    c = Counter()
                    for line in fp:
                        try:
                            # Transform line to dict and cast str numbers to actual numbers
                            read = numeric_cast_dict(keys=col_names, values=line.rstrip().split("\t"))

                            # Filter out ref_id if a select_ref_id list or exclude_ref_id list was provided
                            if select_ref_id and not read["ref_id"] in select_ref_id:
                                raise NanocomporeError("Ref_id not in select list")
                            elif exclude_ref_id and read["ref_id"] in exclude_ref_id:
                                raise NanocomporeError("Ref_id in exclude list")

                            # Filter out reads with high number of invalid kmers if information available
                            if self.__filter_invalid_kmers:
                                if max_invalid_kmers_freq:
                                    invalid_kmers_freq = (read["NNNNN_kmers"]+read["mismatch_kmers"]+read["missing_kmers"])/read["kmers"]
                                    if invalid_kmers_freq > max_invalid_kmers_freq:
                                        raise NanocomporeError("High fraction of invalid kmers ({}%) for read {}".format(round(invalid_kmers_freq*100,2), read["read_id"]))
                                else:
                                    NNNNN_kmers_freq = read["NNNNN_kmers"]/read["kmers"]
                                    max_mismatching_freq = read["mismatch_kmers"]/read["kmers"]
                                    max_missing_freq = read["missing_kmers"]/read["kmers"]
                                    if NNNNN_kmers_freq > max_NNNNN_freq:
                                        raise NanocomporeError("High fraction of NNNNN kmers ({}%) for read {}".format(round(NNNNN_kmers_freq*100,2), read["read_id"]))
                                    elif max_mismatching_freq > max_mismatching_freq:
                                        raise NanocomporeError("High fraction of mismatching kmers ({}%) for read {}".format(round(max_mismatching_freq*100,2), read["read_id"]))
                                    elif max_missing_freq > max_missing_freq:
                                        raise NanocomporeError("High fraction of missing kmers ({}%) for read {}".format(round(max_missing_freq*100,2), read["read_id"]))

                            # Create dict arborescence and save valid reads
                            if not read["ref_id"] in ref_reads:
                                ref_reads[read["ref_id"]] = OrderedDict()
                            if not cond_lab in ref_reads[read["ref_id"]]:
                                ref_reads[read["ref_id"]][cond_lab] = OrderedDict()
                            if not sample_lab in ref_reads[read["ref_id"]][cond_lab]:
                                ref_reads[read["ref_id"]][cond_lab][sample_lab] = []

                            # Fill in list of reads
                            ref_reads[read["ref_id"]][cond_lab][sample_lab].append(read)
                            c ["valid reads"] += 1

                        except NanocomporeError as E:
                            c [str(E)] += 1

                logger.debug("\tCondition:{} Sample:{} {}".format(cond_lab, sample_lab, counter_to_str(c)))
        # Fill in missing condition/sample slots in case
        # a ref_id is missing from one of the eventalign files
        for ref_id in ref_reads.keys():
            for cond_lab, sample_dict in eventalign_fn_dict.items():
                for sample_lab in sample_dict.keys():
                    if not cond_lab in ref_reads[ref_id]:
                        ref_reads[ref_id][cond_lab] = OrderedDict()
                    if not sample_lab in ref_reads[ref_id][cond_lab]:
                        ref_reads[ref_id][cond_lab][sample_lab] = []
        logger.info("\tReferences found in index: {}".format(len(ref_reads)))
        return ref_reads

    def __select_ref(self,
        ref_reads,
        min_coverage,
        min_ref_length,
        downsample_high_coverage):
        """Select ref_id with a minimal coverage in both sample + downsample if needed"""
        valid_ref_reads = OrderedDict()
        c = Counter()
        with Fasta(self._fasta_fn) as fasta:
            for ref_id, ref_dict in ref_reads.items():
                try:
                    # Discard reference transcripts shorter than the threshold
                    assert len(fasta[ref_id]) > min_ref_length
                    valid_dict = OrderedDict()
                    for cond_lab, cond_dict in ref_dict.items():
                        valid_dict[cond_lab] = OrderedDict()
                        for sample_lab, read_list in cond_dict.items():
                            logger.trace(f"Asserting if {ref_id} has enough coverage in {sample_lab}")
                            # Filter out if coverage too low
                            assert len(read_list) >= min_coverage
                            logger.trace(f"ref_id {ref_id} has {len(read_list)} reads in {sample_lab}")
                            # Downsample if coverage too high
                            if downsample_high_coverage and len(read_list) > downsample_high_coverage:
                                read_list = random.sample(read_list, downsample_high_coverage)
                            valid_dict[cond_lab][sample_lab]=read_list


                    # If all valid add to new dict
                    logger.trace(f"ref_id {ref_id} has enough coverage in all samples: keeping it")
                    valid_ref_reads [ref_id] = valid_dict

                    # Save extra info for debug
                    c["valid_ref_id"] += 1
                    for cond_lab, cond_dict in valid_dict.items():
                        for sample_lab, read_list in cond_dict.items():
                            lab = "{} {} Reads".format(cond_lab, sample_lab)
                            c[lab] += len(read_list)

                except AssertionError:
                    logger.trace(f"ref_id {ref_id} does not have enough coverage in at least one sample: discarding it")
                    c["invalid_ref_id"] += 1

        logger.debug(counter_to_str(c))
        logger.info("\tReferences remaining after reference coverage filtering: {}".format(len(valid_ref_reads)))
        return valid_ref_reads
