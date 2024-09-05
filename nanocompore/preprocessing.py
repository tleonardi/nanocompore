import csv
import os
import sqlite3

from collections import Counter
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import closing
from pathlib import Path

import numpy as np

from loguru import logger
from pyfaidx import Fasta

from nanocompore.Transcript import Transcript
from nanocompore.Experiment import Experiment
from nanocompore.Whitelist import Whitelist
from nanocompore.Remora import Remora
from nanocompore.uncalled4 import Uncalled4
from nanocompore.kmer import KmerData
from nanocompore.common import Kit
from nanocompore.common import DROP_KMER_DATA_TABLE_QUERY
from nanocompore.common import DROP_READS_TABLE_QUERY
from nanocompore.common import CREATE_KMER_DATA_TABLE_QUERY
from nanocompore.common import CREATE_READS_TABLE_QUERY
from nanocompore.common import INSERT_KMER_DATA_QUERY
from nanocompore.common import INSERT_READS_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.common import ReadIndexer


class Preprocessor:
    """
    The preprocessor takes data coming from
    the chosen resquiggler (described in the
    configuration file used for running the
    experiment) and prepares it in a common
    format in a SQLite DB to be used in the
    subsequest analysis step.
    """
    def __init__(self, config):
        logger.info("Initializing the preprocessor")
        self._config = config
        self._experiment = Experiment(config)
        self._whitelist_transcripts()
        db = self._config.get_kmer_data_db()
        self._db_conn = sqlite3.connect(db)
        self._db_cursor = self._db_conn.cursor()
        self._db_cursor.execute(DROP_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(CREATE_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(DROP_READS_TABLE_QUERY)
        self._db_cursor.execute(CREATE_READS_TABLE_QUERY)
        self._read_indexer = ReadIndexer()
        self._sample_ids = self._experiment.get_sample_ids()
        # The preprocessor will use the main process
        # and will spawn n - 1 worker processes in
        # order to make sure we use the number of
        # processes specified by the user in the config.
        self._worker_processes = self._config.get_nthreads() - 1


    def _whitelist_transcripts(self):
        logger.info("Starting to whitelist the reference IDs")
        whitelist = Whitelist(self._experiment, self._config)

        self._valid_transcripts = whitelist.ref_id_list


    def _write_to_db(self, ref_id, kmer_data_generator):
        self._db_cursor.executemany(
                INSERT_KMER_DATA_QUERY,
                list(self._process_rows_for_writing(ref_id, kmer_data_generator)))
        self._db_conn.commit()


    def _process_rows_for_writing(self, ref_id, kmer_data_generator):
        for kmer_data in kmer_data_generator:
            new_mappings = self._read_indexer.add_reads(kmer_data.reads)
            self._db_cursor.executemany(INSERT_READS_QUERY, new_mappings)

            read_ids = self._read_indexer.get_ids(kmer_data.reads)
            read_ids = np.array(read_ids, dtype=READ_ID_TYPE)

            sample_ids = [self._sample_ids[label]
                          for label in kmer_data.sample_labels]
            sample_ids = np.array(sample_ids, dtype=SAMPLE_ID_TYPE)

            yield (ref_id,
                   kmer_data.pos,
                   kmer_data.kmer,
                   sample_ids.tobytes(),
                   read_ids.tobytes(),
                   kmer_data.intensity.tobytes(),
                   kmer_data.sd.tobytes(),
                   kmer_data.dwell.tobytes())


    def __del__(self):
        if '_db_cursor' in self.__dict__:
            self._db_cursor.close()
        if '_db_conn' in self.__dict__:
            self._db_conn.close()


    # Note: when we use the preprocessor's instance
    # methods with multiprocessing the library will
    # try to pickle the instance object so that it
    # can send the method's self to the worker process.
    # Since the database connection cannot be pickled
    # and we don't want to use it in the worker processes,
    # but only in the main process, we override the
    # getstate/setstate methods to exclude the db
    # connection and other potentially big variables
    # from the set of variables that will be sent to
    # the worker processes.
    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['_db_conn']
        del state['_db_cursor']
        del state['_read_indexer']
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.__dict__.update(state)


class RemoraPreprocessor(Preprocessor):
    """
    Uses ONT's Remora internally to resquiggle
    the input and write the results to an SQLite
    database for later analysis.
    """

    def __init__(self, config):
        super().__init__(config)


    def __call__(self):
        with ProcessPoolExecutor(max_workers=self._worker_processes) as executor:
            futures = [executor.submit(self._resquiggle, ref_id)
                       for ref_id in self._valid_transcripts]
            for future in as_completed(futures):
                ref_id, kmer_data_generator = future.result()
                self._write_to_db(ref_id, kmer_data_generator)


    def _resquiggle(self, ref_id):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])

        remora = Remora(self._experiment,
                        self._config,
                        ref_id=ref_id,
                        start=0,
                        end=len(ref_seq),
                        seq=ref_seq,
                        strand='+')

        return ref_id, list(remora.kmer_data_generator())


class Uncalled4Preprocessor(Preprocessor):
    """
    Preprocess a bam file produced by the
    uncalled4 resquiggler's align command.
    The signal data would be transfered
    to an SQLite DB to be used by Nanocompore
    for later analysis.
    """

    def __init__(self, config):
        super().__init__(config)


    def __call__(self):
        with ProcessPoolExecutor(max_workers=self._worker_processes) as executor:
            futures = [executor.submit(self._resquiggle, ref_id)
                       for ref_id in self._valid_transcripts]
            for future in as_completed(futures):
                ref_id, kmer_data_generator = future.result()
                self._write_to_db(ref_id, kmer_data_generator)


    def _resquiggle(self, ref_id):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])
        uncalled4 = Uncalled4(self._experiment,
                              self._config,
                              ref_id,
                              ref_seq)
        return ref_id, list(uncalled4.kmer_data_generator())

