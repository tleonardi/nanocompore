import csv
import os
import sqlite3
import json

import multiprocessing as mp

from collections import Counter
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import closing
from pathlib import Path

import pysam
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
from nanocompore.common import DROP_METADATA_TABLE_QUERY
from nanocompore.common import DROP_KMER_DATA_TABLE_QUERY
from nanocompore.common import DROP_READS_TABLE_QUERY
from nanocompore.common import DROP_READS_ID_INDEX_QUERY
from nanocompore.common import DROP_KMERS_INDEX_QUERY
from nanocompore.common import DROP_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import CREATE_METADATA_TABLE_QUERY
from nanocompore.common import CREATE_KMER_DATA_TABLE_QUERY
from nanocompore.common import CREATE_READS_TABLE_QUERY
from nanocompore.common import CREATE_READS_ID_INDEX_QUERY
from nanocompore.common import CREATE_KMERS_INDEX_QUERY
from nanocompore.common import CREATE_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import INSERT_KMER_DATA_QUERY
from nanocompore.common import INSERT_READS_QUERY
from nanocompore.common import INSERT_TRANSCRIPTS_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.common import Indexer
from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.common import get_measurement_type
from nanocompore.common import DB_METADATA_RESQUIGGLER_KEY
from nanocompore.common import DB_METADATA_READ_ID_TYPE_KEY
from nanocompore.common import DB_METADATA_SAMPLE_ID_TYPE_KEY
from nanocompore.common import DB_METADATA_MEASUREMENT_TYPE_KEY
from nanocompore.common import DB_METADATA_SAMPLE_LABELS_KEY
from nanocompore.common import DB_METADATA_CONDITION_SAMPLES_KEY


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

        # self._db_conn = sqlite3.connect(db, isolation_level=None)
        # self._db_conn.execute('PRAGMA synchronous = OFF')
        # self._db_conn.execute('PRAGMA journal_mode = OFF')
        # # self._db_conn.execute('PRAGMA cache_size = -200000')
        self._db_conn = self._get_db()
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
        with closing(self._db_conn.cursor()) as cursor:
            cursor.execute(DROP_KMER_DATA_TABLE_QUERY)
            cursor.execute(CREATE_KMER_DATA_TABLE_QUERY)
            cursor.execute(DROP_READS_TABLE_QUERY)
            cursor.execute(CREATE_READS_TABLE_QUERY)
            cursor.execute(DROP_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(DROP_METADATA_TABLE_QUERY)
            cursor.execute(CREATE_METADATA_TABLE_QUERY)
            cursor.execute(DROP_READS_ID_INDEX_QUERY)
            cursor.execute(DROP_KMERS_INDEX_QUERY)
        self._write_metadata()


    def _create_db_indices(self):
        with closing(self._db_conn.cursor()) as cursor:
            cursor.execute(CREATE_READS_ID_INDEX_QUERY)
            cursor.execute(CREATE_KMERS_INDEX_QUERY)


    def _write_metadata(self):
        resquiggler = self._config.get_resquiggler()
        condition_samples = {cond: self._experiment.condition_to_samples(cond)
                             for cond in self._experiment.get_condition_labels()}
        metadata = {
            DB_METADATA_RESQUIGGLER_KEY: resquiggler,
            DB_METADATA_READ_ID_TYPE_KEY: READ_ID_TYPE.__name__,
            DB_METADATA_SAMPLE_ID_TYPE_KEY: SAMPLE_ID_TYPE.__name__,
            DB_METADATA_MEASUREMENT_TYPE_KEY: get_measurement_type(resquiggler).__name__,
            DB_METADATA_SAMPLE_LABELS_KEY: ','.join(self._experiment.get_sample_labels()),
            DB_METADATA_CONDITION_SAMPLES_KEY: json.dumps(condition_samples),
        }
        with closing(self._db_conn.cursor()) as cursor:
            for key, value in metadata.items():
                cursor.execute("INSERT INTO metadata (key, value) VALUES (?, ?)", (key, str(value)))


    def _get_references_from_bams(self):
        logger.info("Getting references from the BAMs.")
        references = set()
        for condition_def in self._config.get_data().values():
            for sample, sample_def in condition_def.items():
                bam = pysam.AlignmentFile(sample_def['bam'], "rb")
                references.update(bam.references)
        logger.info(f"Found {len(references)} references.")
        return references


    def _get_db(self, db=None):
        if not db:
            db = self._config.get_kmer_data_db()
        conn = sqlite3.connect(db, isolation_level=None)
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = OFF')
        # Use largest possible page size of 65kb
        conn.execute('PRAGMA page_size = 65536')
        # Use cache of 15000xpage_size. That's about
        # 952mb of cache.
        conn.execute('PRAGMA cache_size = 2000')
        # conn.execute('PRAGMA cache_size = -400000')
        return conn


    def _write_kmer_rows(self, connection, rows):
        with closing(connection.cursor()) as cursor:
            cursor.execute("begin")
            cursor.executemany(INSERT_KMER_DATA_QUERY, rows)
            cursor.execute("commit")
        # self._db_conn.commit()


    def _process_rows_for_writing(self, ref_id, kmers_data, read_invalid_ratios):
        with closing(self._db_conn.cursor()) as cursor:
            cursor.execute(INSERT_TRANSCRIPTS_QUERY,
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
        reads_data = [(read, idx, read_invalid_ratios[read])
                      for read, idx in new_mappings]
        with closing(self._db_conn.cursor()) as cursor:
            cursor.execute("begin")
            cursor.executemany(INSERT_READS_QUERY, reads_data)
            cursor.execute("commit")

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
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.__dict__.update(state)


    def _get_reads_invalid_kmer_ratio(self, kmers):
        """
        Calculate the ratio of missing kmers
        in the read. It takes the kmers with
        min and max position to determine the
        read length.

        This is the method employed for Uncalled4 and Remora.
        For eventalign there's a custom method, that takes
        in consideration the richer information provided
        by the resquiggler.
        """
        read_counts = defaultdict(lambda: 0)
        read_ends = defaultdict(lambda: (np.inf, -1))

        for kmer in kmers:
            for read in kmer.reads:
                read_counts[read] += 1
                curr_range = read_ends[read]
                start = min(kmer.pos, curr_range[0])
                end = max(kmer.pos, curr_range[1])
                read_ends[read] = (start, end)
        return {read: self._calc_invalid_ratio(read_ends[read], count)
                for read, count in read_counts.items()}


    def _calc_invalid_ratio(self, ends, valid):
        length = ends[1] - ends[0] + 1
        return (length - valid)/length


class RemoraPreprocessor(Preprocessor):
    """
    Uses ONT's Remora internally to resquiggle
    the input and write the results to an SQLite
    database for later analysis.
    """

    def __init__(self, config):
        super().__init__(config)
        self._references = self._get_references_from_bams()


    def __call__(self):
        with ProcessPoolExecutor(max_workers=self._worker_processes) as executor:
            futures = [executor.submit(self._resquiggle, ref_id)
                       for ref_id in self._references]
            for future in as_completed(futures):
                ref_id, kmers, read_invalid_ratios = future.result()
                rows = self._process_rows_for_writing(ref_id, kmers, read_invalid_ratios)
                with closing(self._get_db()) as conn:
                    self._write_kmer_rows(conn, rows)
        self._create_db_indices()


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
        kmers = list(remora.kmer_data_generator())
        read_invalid_ratios = self._get_reads_invalid_kmer_ratio(kmers)
        return ref_id, kmers, read_invalid_ratios


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
        self._references = self._get_references_from_bams()


    def __call__(self):
        with ProcessPoolExecutor(max_workers=self._worker_processes) as executor:
            futures = [executor.submit(self._resquiggle, ref_id)
                       for ref_id in self._references]
            for future in as_completed(futures):
                ref_id, kmers, read_invalid_ratios = future.result()
                rows = self._process_rows_for_writing(ref_id, kmers, read_invalid_ratios)
                with closing(self._get_db()) as conn:
                    self._write_kmer_rows(conn, rows)
        self._create_db_indices()


    def _resquiggle(self, ref_id):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])
        uncalled4 = Uncalled4(self._experiment,
                              self._config,
                              ref_id,
                              ref_seq)
        kmers = list(uncalled4.kmer_data_generator())
        read_invalid_ratios = self._get_reads_invalid_kmer_ratio(kmers)
        return ref_id, kmers, read_invalid_ratios


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

        logger.info("Reviewing input collapsed dbs and transferring read id mappings to the new db.")

        sample_defs = {sample: sample_def
                       for condition_def in self._config.get_data().values()
                       for sample, sample_def in condition_def.items()}

        references = {}
        total_reads = 0
        sample_read_offsets = {}
        current_offset = 0
        for sample, sample_def in sample_defs.items():
            sample_read_offsets[sample] = current_offset

            db_in_path = sample_def['eventalign_db']
            with closing(sqlite3.connect(db_in_path)) as conn,\
                 closing(conn.cursor()) as db_in,\
                 closing(self._db_conn.cursor()) as db_out:

                for row in db_in.execute("SELECT reference FROM transcripts").fetchall():
                    references[row[0]] = True

                # Get all reads from the collapsed sample db,
                # offset the ids (so that they don't overlap with those
                # of other samples) and move the offsetted mappings to
                # the merged db.
                reads = db_in.execute("SELECT * FROM reads").fetchall()
                total_reads += len(reads)
                offsetted_reads = [(row[0], row[1] + current_offset, row[2]) for row in reads]
                db_out.execute("begin")
                db_out.executemany(INSERT_READS_QUERY, offsetted_reads)
                db_out.execute("commit")

                max_read_id = db_in.execute("SELECT MAX(id) FROM reads").fetchone()[0]
                current_offset += max_read_id

        references = list(references.keys())
        logger.info(f"{len(references)} references found.")
        logger.info(f"{total_reads} read mappings were transfered.")

        logger.info(f"Starting to read kmer data from sample dbs and merging it.")

        # Merge the information from the collapsed eventalign_dbs

        # We start worker processes and put all
        # references in a task queue, from where
        # the workers will pick them up and
        # process them. When a worker finishes
        # the processing of a reference, it will
        # acquire a lock and save the results to
        # the database. We use this more cumbersome
        # approach with manually creating proceses,
        # queues and locks instead of using a
        # ProcessPoolExecutor, because for references
        # with high coverage the writing can be slow
        # and the workers can outpace the main consumer
        # process, which leads accumulation of big
        # results from the pool and memory spikes.

        manager = mp.Manager()
        lock = manager.Lock()
        task_queue = manager.Queue()
        current_transcript_id = manager.Value('i', 1)
        processed_transcripts = manager.Value('i', 0)

        workers = [mp.Process(target=self._read_and_merge_samples,
                              args=(i,
                                    task_queue,
                                    lock,
                                    sample_read_offsets,
                                    current_transcript_id,
                                    processed_transcripts,
                                    len(references)))
                   for i in range(self._worker_processes)]

        for worker in workers:
            worker.start()

        for ref_id in references:
            task_queue.put(ref_id)

        # add poison pills to kill the workers
        for _ in range(self._worker_processes):
            task_queue.put(None)

        for worker in workers:
            worker.join()

        logger.info("Merging tmp databases.")

        # Merge tmp databases
        with closing(self._get_db()) as conn,\
             closing(conn.cursor()) as cursor:
            for i in range(self._worker_processes):
                tmp_db = self._config.get_kmer_data_db() + f".{i}"
                cursor.execute(f"ATTACH '{tmp_db}' as tmp")
                cursor.execute("BEGIN")
                cursor.execute("INSERT INTO transcripts SELECT * FROM tmp.transcripts")
                cursor.execute("INSERT INTO kmer_data SELECT * FROM tmp.kmer_data")
                conn.commit()
                cursor.execute("DETACH DATABASE tmp")
                # Delete the tmp db
                Path(tmp_db).unlink()

        logger.info("Creating database indices.")
        self._create_db_indices()


    def _read_and_merge_samples(self,
                                idx,
                                task_queue,
                                lock,
                                sample_read_offsets,
                                current_transcript_id,
                                processed_transcripts,
                                num_transcripts):
        sample_dbs = {sample: sample_def['eventalign_db']
                      for condition_def in self._config.get_data().values()
                      for sample, sample_def in condition_def.items()}

        tmp_db_out = self._config.get_kmer_data_db() + f".{idx}"
        with closing(self._get_db(tmp_db_out)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(DROP_KMER_DATA_TABLE_QUERY)
            cursor.execute(CREATE_KMER_DATA_TABLE_QUERY)
            cursor.execute(DROP_READS_TABLE_QUERY)
            cursor.execute(CREATE_READS_TABLE_QUERY)
            cursor.execute(DROP_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)
            conn.commit()

        while True:
            ref_id = task_queue.get()
            if ref_id == None:
                break

            pos_data = defaultdict(list)
            for sample, db in sample_dbs.items():
                with closing(sqlite3.connect(db)) as conn,\
                     closing(conn.cursor()) as cursor:
                    transcript_id = cursor.execute("SELECT id FROM transcripts WHERE reference = ?",
                                                   (ref_id,)).fetchone()[0]
                    rows = cursor.execute("SELECT * FROM kmer_data WHERE transcript_id = ?",
                                          (transcript_id,)).fetchall()
                    for row in rows:
                        pos = row[1]
                        pos_data[pos].append((sample, row))

            kmers = self._merge_kmer_data_from_samples(pos_data, sample_read_offsets)
            if not kmers:
                continue

            rows = self._process_rows_for_writing(ref_id, kmers, None)

            lock.acquire()
            transcript_id = current_transcript_id.value
            current_transcript_id.value += 1
            lock.release()

            with closing(self._get_db(tmp_db_out)) as conn,\
                 closing(conn.cursor()) as cursor:
                cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                               (ref_id, transcript_id))
                self._write_kmer_rows(conn, rows)
            lock.acquire()
            processed_transcripts.value += 1
            logger.info(f"Data for reference {processed_transcripts.value}/{num_transcripts} {ref_id} has been saved to the tmp database.")
            lock.release()


    def _merge_kmer_data_from_samples(self, pos_data, sample_read_offsets):
        kmers = []
        for pos, rows in pos_data.items():
            kmer = rows[0][1][2]
            samples = np.concatenate([np.repeat(sample, len(row[3]) / np.dtype(READ_ID_TYPE).itemsize)
                                      for sample, row in rows])
            read_ids = [idx + sample_read_offsets[sample]
                        for idx, sample in zip(self._merge_col_from_samples(3, rows, READ_ID_TYPE),
                                               samples)]
            read_ids = np.array(read_ids, dtype=READ_ID_TYPE)

            intensity = self._merge_col_from_samples(4, rows, EVENTALIGN_MEASUREMENT_TYPE)
            sd = self._merge_col_from_samples(5, rows, EVENTALIGN_MEASUREMENT_TYPE)
            dwell = self._merge_col_from_samples(6, rows, EVENTALIGN_MEASUREMENT_TYPE)

            kmers.append(KmerData(pos, kmer, samples, read_ids, intensity, sd, dwell, None, self._experiment))
        return kmers


    def _merge_col_from_samples(self, col, rows, dtype):
        joined_bytes = bytes().join([row[col] for _, row in rows])
        return np.frombuffer(joined_bytes, dtype=dtype)


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
                                        self._get_intermediary_db_name(sample),
                                        self._config.get_kit()))
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


    # Since the kmer data we read from the eventalign collapsed db
    # already has internal integer indices instead of uuids,
    # we need a different processing that will not try to map them again.
    # Also, in this case the read mappings are already transfered from
    # the eventalign collapsed db and we don't need to write them here.
    def _process_rows_for_writing(self, ref_id, kmers_data, read_invalid_ratios=None):
        processed_kmers = []
        for kmer_data in kmers_data:
            sample_ids = [self._sample_ids[label]
                          for label in kmer_data.sample_labels]
            sample_ids = np.array(sample_ids, dtype=SAMPLE_ID_TYPE)

            proc_kmer = (self._current_transcript_id,
                         kmer_data.pos,
                         kmer_data.kmer,
                         sample_ids.tobytes(),
                         kmer_data.reads.tobytes(),
                         kmer_data.intensity.tobytes(),
                         kmer_data.sd.tobytes(),
                         kmer_data.dwell.tobytes())
            processed_kmers.append(proc_kmer)
        return processed_kmers


def collapse_eventalign(params):
    sample, eventalign, fasta_ref, output, kit = params
    EventalignCollapser(eventalign, fasta_ref, output, kit)()
    return sample, output

