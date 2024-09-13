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
from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.kmer import KmerData
from nanocompore.common import Kit
from nanocompore.common import DROP_KMER_DATA_TABLE_QUERY
from nanocompore.common import DROP_READS_TABLE_QUERY
from nanocompore.common import DROP_READS_ID_INDEX_QUERY
from nanocompore.common import DROP_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import CREATE_KMER_DATA_TABLE_QUERY
from nanocompore.common import CREATE_READS_TABLE_QUERY
from nanocompore.common import CREATE_READS_ID_INDEX_QUERY
from nanocompore.common import CREATE_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import INSERT_KMER_DATA_QUERY
from nanocompore.common import INSERT_READS_QUERY
from nanocompore.common import INSERT_TRANSCRIPTS_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.common import Indexer
from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.common import get_reads_invalid_kmer_ratio


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
        self._setup_db()

        self._sample_ids = self._experiment.get_sample_ids()

        self._current_transcript_id = 1
        self._current_read_id = 1

        # The preprocessor will use the main process
        # and will spawn n - 1 worker processes in
        # order to make sure we use the number of
        # processes specified by the user in the config.
        self._worker_processes = self._config.get_nthreads() - 1


    def _setup_db(self):
        self._db_cursor.execute(DROP_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(CREATE_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(DROP_READS_TABLE_QUERY)
        self._db_cursor.execute(CREATE_READS_TABLE_QUERY)
        self._db_cursor.execute(DROP_READS_ID_INDEX_QUERY)
        self._db_cursor.execute(CREATE_READS_ID_INDEX_QUERY)
        self._db_cursor.execute(DROP_TRANSCRIPTS_TABLE_QUERY)
        self._db_cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)


    def _whitelist_transcripts(self):
        logger.info("Starting to whitelist the reference IDs")
        whitelist = Whitelist(self._experiment, self._config)

        self._valid_transcripts = whitelist.ref_id_list


    def _write_to_db(self, ref_id, kmers_data):
        rows = list(self._process_rows_for_writing(ref_id, kmers_data))
        self._db_cursor.executemany(INSERT_KMER_DATA_QUERY, rows)
        self._db_conn.commit()


    def _process_rows_for_writing(self, ref_id, kmers_data):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_len = len(fasta_fh[ref_id])

        read_invalid_kmer_ratios = get_reads_invalid_kmer_ratio(kmers_data, ref_len)

        self._db_cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                                (ref_id, self._current_transcript_id))

        # We assume that all reads for a single
        # transcript would be processed in a single
        # call, so we don't need to store all
        # read ids for all transcripts for the
        # whole execution of the preprocessor.
        # Instead we create a new indexer for
        # each transcript, but we make sure
        # that ids are not repeated by making
        # the new indexer starting from the last
        # id of the previous indexer.
        read_indexer = Indexer(initial_index=self._current_read_id)
        all_reads = {read
                     for kmer in kmers_data
                     for read in kmer.reads}
        new_mappings = read_indexer.add(all_reads)
        reads_data = [(read, idx, read_invalid_kmer_ratios[read])
                      for read, idx in new_mappings]
        self._db_cursor.executemany(INSERT_READS_QUERY, reads_data)

        for kmer_data in kmers_data:
            read_ids = read_indexer.get_ids(kmer_data.reads)
            read_ids = np.array(read_ids, dtype=READ_ID_TYPE)

            sample_ids = [self._sample_ids[label]
                          for label in kmer_data.sample_labels]
            sample_ids = np.array(sample_ids, dtype=SAMPLE_ID_TYPE)

            yield (self._current_transcript_id,
                   kmer_data.pos,
                   kmer_data.kmer,
                   sample_ids.tobytes(),
                   read_ids.tobytes(),
                   kmer_data.intensity.tobytes(),
                   kmer_data.sd.tobytes(),
                   kmer_data.dwell.tobytes())

        self._current_transcript_id += 1
        self._current_read_id = read_indexer.current_id


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
                ref_id, kmers_data = future.result()
                self._write_to_db(ref_id, kmers_data)


    def _resquiggle(self, ref_id):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])
        uncalled4 = Uncalled4(self._experiment,
                              self._config,
                              ref_id,
                              ref_seq)
        return ref_id, list(uncalled4.kmer_data_generator())


class EventalignPreprocessor(Preprocessor):
    """
    Takes the output of Nanopolish or F5C eventalign
    command and prepares the data for nanocompore
    by collapsing it and storing it in an SQLite DB
    for later analysis.
    """

    def __init__(self, config):
        super().__init__(config)

        self._validate_eventalign_input()


    def __call__(self):
        self._reuse_collapsed_files()

        # If any of the samples has a raw eventalign
        # file as input, collapse it.
        if any('eventalign_db' not in sample_def
               for condition_def in self._config.get_data().values()
               for sample_def in condition_def.values()):
            self._collapse_eventaligns()

        sample_dbs = {sample: sqlite3.connect(sample_def['eventalign_db'])
                      for condition_def in self._config.get_data().values()
                      for sample, sample_def in condition_def.items()}

        sample_cursors = {sample: db.cursor()
                          for sample, db in sample_dbs.items()}

        sample_read_indices = {}
        for sample, cursor in sample_cursors.items():
            res = cursor.execute("SELECT * FROM reads")
            sample_read_indices[sample] = {row[1]: row[0]
                                           for row in res.fetchall()}

        # Merge the information from the collapsed eventalign_dbs
        for ref_id in self._valid_transcripts:
            pos_data = defaultdict(list)
            for sample, cursor in sample_cursors.items():
                transcript_id = cursor.execute("SELECT id FROM transcripts WHERE reference = ?", (ref_id,)).fetchone()[0]
                rows = cursor.execute("SELECT * FROM kmer_data WHERE transcript_id = ?", (transcript_id,)).fetchall()
                for row in rows:
                    pos = row[1]
                    pos_data[pos].append((sample, row))
            kmers_list = self._merge_kmer_data_from_samples(pos_data, sample_read_indices)
            self._write_to_db(ref_id, kmers_list)


    def _merge_kmer_data_from_samples(self, pos_data, sample_read_indices):
        kmers = []
        for pos, rows in pos_data.items():
            kmer = rows[0][1][2]
            samples = np.concatenate([np.repeat(sample, len(np.frombuffer(row[3], dtype=READ_ID_TYPE)))
                                      for sample, row in rows])
            read_ids = self._merge_reads_from_rows(3, rows, sample_read_indices)
            intensity = self._merge_measurements_from_rows(4, rows, EVENTALIGN_MEASUREMENT_TYPE)
            sd = self._merge_measurements_from_rows(5, rows, EVENTALIGN_MEASUREMENT_TYPE)
            dwell = self._merge_measurements_from_rows(6, rows, EVENTALIGN_MEASUREMENT_TYPE)

            kmers.append(KmerData(pos, kmer, samples, read_ids, intensity, sd, dwell, self._experiment))
        return kmers


    def _merge_reads_from_rows(self, col, rows, sample_read_indices):
        return np.array([sample_read_indices[sample][index]
                         for sample, row in rows
                         for index in np.frombuffer(row[col], dtype=READ_ID_TYPE)])


    def _merge_measurements_from_rows(self, col, rows, dtype):
        return np.array([value
                         for sample, row in rows
                         for value in np.frombuffer(row[col], dtype=dtype)])


    def _reuse_collapsed_files(self):
        for condition_def in self._config.get_data().values():
            for sample, sample_def in condition_def.items():
                if 'eventalign_db' in sample_def:
                    continue

                intermediary_db = self._get_intermediary_db_name(sample)
                if os.path.isfile(intermediary_db):
                    logger.warning(f"Sample {sample} seems to already have been collapsed. " +\
                                   f"Reusing the existing file: {intermediary_db}")
                    sample_def['eventalign_db'] = intermediary_db


    def _collapse_eventaligns(self):
        noncollapsed = {sample: sample_def['eventalign_tsv']
                        for condition_def in self._config.get_data().values()
                        for sample, sample_def in condition_def.items()
                        if 'eventalign_db' not in sample_def}
        logger.info(f"{len(noncollapsed)} samples have input eventalign files that need to be collapsed. Collapsing them now.")

        num_processes = min(self._worker_processes, len(noncollapsed))
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            futures = [executor.submit(collapse_eventalign,
                                       (sample,
                                        file,
                                        self._config.get_fasta_ref(),
                                        self._get_intermediary_db_name(sample)))
                       for sample, file in noncollapsed.items()]
            for future in as_completed(futures):
                sample, db = future.result()
                logger.info(f"Input eventalign for sample {sample} has been collapsed and saved at {db}")
                condition = self._experiment.sample_to_condition(sample)
                self._config.get_data()[condition][sample]['eventalign_db'] = db


    def _get_intermediary_db_name(self, sample):
        kmer_db = self._config.get_kmer_data_db()
        path = Path(kmer_db)
        return str(path.with_name(path.stem + '_' + sample + path.suffix))


    def _validate_eventalign_input(self):
        for condition in self._config.get_data().values():
            for sample in condition.values():
                if 'eventalign_tsv' not in sample and 'eventalign_db' not in sample:
                    raise ValueError('When using the "eventalign" preprocessor each sample must contain either the field "eventalign_tsv" with a path to the eventalign tsv file or "eventalign_db" with the alreday collapsed eventalign data.')


def collapse_eventalign(params):
    sample, eventalign, fasta_ref, output = params
    EventalignCollapser(eventalign, fasta_ref, output)()
    return sample, output

