# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import logging
from loguru import logger
import random
import sqlite3

# Third party
import numpy as np
from tqdm import tqdm
from pyfaidx import Fasta

# Local package
from nanocompore.common import *
from nanocompore.DatabaseWrapper import DatabaseWrapper

# Set global random seed
downsample_random_seed = 42
random.seed(downsample_random_seed)

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist(object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
                 db_path,
                 sample_dict,
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
        Generate a whitelist of reads that fulfill filtering criteria
        Args:
        * db_path
            Path to the SQLite database file with event-aligned read/kmer data
        * sample_dict
            Dictionary containing lists of (unique) sample names, grouped by condition
            example d = {"control": ["C1", "C2"], "treatment": ["T1", "T2"]}
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

        check_sample_dict(sample_dict)

        # Create look-up dict. of conditions
        cond_dict = {}
        for cond, samples in sample_dict.items():
            for sample in samples:
                cond_dict[sample] = cond

        # Test if Fasta can be opened
        try:
            with Fasta(fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        # Database interaction
        with DatabaseWrapper(db_path) as db:
            db_samples = db.get_samples(sample_dict)

            # Set up filters by adding conditions for DB query
            select = ["reads.id AS readid", "sampleid", "transcriptid", "transcripts.name AS transcriptname"]
            where = []
            # Get reads only from a subset of samples?
            if len(cond_dict) != len(db_samples):
                where = ["sampleid IN (%s)" % ", ".join(map(str, db_samples))]

            if select_ref_id:
                select.append("reads.name AS readname")
                where.append("readname IN ('%s')" % "', '".join(select_ref_id))
            elif exclude_ref_id:
                select.append("reads.name AS readname")
                where.append("readname NOT IN ('%s')" % "', '".join(exclude_ref_id))

            if max_invalid_kmers_freq is not None:
                if max_invalid_kmers_freq < 1.0:
                    select.append("1.0 - CAST(valid_kmers AS REAL) / kmers AS invalid_freq")
                    where.append(f"invalid_freq <= {max_invalid_kmers_freq}")
            else:
                if max_NNNNN_freq < 1.0:
                    select.append("CAST(NNNNN_kmers AS REAL) / kmers AS NNNNN_freq")
                    where.append(f"NNNNN_freq <= {max_NNNNN_freq}")
                if max_mismatching_freq < 1.0:
                    select.append("CAST(mismatch_kmers AS REAL) / kmers AS mismatch_freq")
                    where.append(f"mismatch_freq <= {max_mismatching_freq}")
                if max_missing_freq < 1.0:
                    select.append("CAST(missing_kmers AS REAL) / kmers AS missing_freq")
                    where.append(f"missing_freq <= {max_missing_freq}")

            query = "SELECT %s FROM reads LEFT JOIN transcripts ON transcriptid = transcripts.id" % \
                ", ".join(select)
            if where:
                query += " WHERE %s" % " AND ".join(where)

            # dict. structure: transcript -> condition -> sample -> list of reads
            ref_reads = {}
            logger.info("Querying reads from DB")
            try:
                db.cursor.execute(query)
                for row in db.cursor:
                    read_id = row["readid"]
                    sample_id = row["sampleid"]
                    condition = cond_dict[db_samples[sample_id]]
                    ref_id = row["transcriptname"]
                    ref_reads.setdefault(ref_id, {}).setdefault(condition, {}).\
                        setdefault(sample_id, []).append(read_id)
            except Exception:
                logger.error("Error querying reads from DB")
                raise Exception

        # Filtering at transcript level
        logger.info("Filtering out references with low coverage")
        self.ref_reads = self.__select_ref(
            ref_reads=ref_reads,
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
