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

# Local package
from nanocompore.common import *
from nanocompore.DataStore import DataStore_transcript

# Set global random seed
downsample_random_seed = 42
random.seed(downsample_random_seed)

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist(object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
                 sample_dict,
                 min_coverage = 10,
                 max_invalid_kmers_freq = 0.1,
                 max_NNNNN_freq = 0.1,
                 max_mismatching_freq = 0.1,
                 max_missing_freq = 0.1):
        # TODO: add 'min_coverage' equivalent for conditions, not samples
        """
        Generate a whitelist of reads that fulfill filtering criteria
        Args:
        * sample_dict
            Dictionary containing lists of sample ids and names, grouped by condition
            example d = {"control": [(1, "C1"), (2, "C2")], "treatment": [(3, "T1"), (4, "T2")]}
        * min_coverage
            minimal coverage required in each sample
        * max_invalid_kmers_freq
            maximum frequency of NNNNN, mismatching and missing kmers in reads
            If None, then 'max_NNNNN_freq', 'max_mismatching_freq', 'max_missing_freq' will be used instead
        * max_NNNNN_freq
            maximum frequency of NNNNN kmers in reads (1 to deactivate)
        * max_mismatching_freq
            maximum frequency of mismatching kmers in reads (1 to deactivate)
        * max_missing_freq
            maximum frequency of missing kmers in reads (1 to deactivate)
        """
        # Total number of samples:
        self._n_samples = sum(map(len, sample_dict.values()))

        where = None
        if max_invalid_kmers_freq is not None:
            if max_invalid_kmers_freq < 1.0:
                where = f"1.0 - CAST(valid_kmers AS REAL) / kmers <= {max_invalid_kmers_freq}"
        else:
            conditions = []
            if max_NNNNN_freq < 1.0:
                conditions.append(f"CAST(NNNNN_kmers AS REAL) / kmers <= {max_NNNNN_freq}")
            if max_mismatching_freq < 1.0:
                conditions.append(f"CAST(mismatch_kmers AS REAL) / kmers <= {max_mismatching_freq}")
            if max_missing_freq < 1.0:
                conditions.append(f"CAST(missing_kmers AS REAL) / kmers <= {max_missing_freq}")
            if conditions:
                where = " AND ".join(conditions)

        self._update_query = "UPDATE reads SET pass_filter = 1"
        if where:
            self._update_query += " WHERE " + where

        self._min_coverage = min_coverage


    def __call__(self, db_path):
        """Filter reads in a transcript database, return whether the transcript meets criteria"""
        with DataStore_transcript(db_path, "", 0) as db:
            # filter reads:
            try:
                db.cursor.execute("ALTER TABLE reads ADD COLUMN pass_filter INTEGER DEFAULT 0")
            except sqlite3.OperationalError as error:
                if error.args[0] == "duplicate column name: pass_filter":
                    # column exists from previous round of filtering - reset values to 0:
                    db.cursor.execute("UPDATE reads SET pass_filter = 0")
                else:
                    raise
            db.cursor.execute(self._update_query)
            # count reads that passed the filter:
            if self._min_coverage:
                subquery = "SELECT sampleid, COUNT(*) AS n_reads FROM reads WHERE pass_filter = 1 GROUP BY sampleid"
                sql = f"SELECT COUNT(*) AS n_samples, MIN(n_reads) AS min_reads FROM ({subquery})"
                db.cursor.execute(sql)
                row = db.cursor.fetchone()
                if (row["n_samples"] < self._n_samples) or (row["min_reads"] < self._min_coverage):
                    return False
        return True
