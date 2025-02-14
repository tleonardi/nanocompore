import csv
import os
import sqlite3

import multiprocessing as mp

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from contextlib import closing
from pathlib import Path

import pysam
import numpy as np

from loguru import logger
from pyfaidx import Fasta

from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.common import Indexer
from nanocompore.common import Kit
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.database import INSERT_READS_QUERY
from nanocompore.database import INSERT_TRANSCRIPTS_QUERY
from nanocompore.database import PreprocessingDB
from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.kmer import KmerData
from nanocompore.remora_wrapper import Remora
from nanocompore.uncalled4 import Uncalled4


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

        self._db = PreprocessingDB(config.get_preprocessing_db(), config)
        self._db.setup()

        self._sample_ids = self._config.get_sample_ids()

        self._current_transcript_id = 1
        self._current_read_id = 1

        # The preprocessor will use the main process
        # and will spawn n - 1 worker processes in
        # order to make sure we use the number of
        # processes specified by the user in the config.
        self._worker_processes = self._config.get_nthreads() - 1


    def _get_references_from_bams(self):
        logger.info("Getting references from the BAMs.")
        references = set()
        for condition_def in self._config.get_data().values():
            for sample, sample_def in condition_def.items():
                bam = pysam.AlignmentFile(sample_def['bam'], "rb")
                for line in bam.get_index_statistics():
                    if line.mapped > 0:
                        references.add(line.contig)
        logger.info(f"Found {len(references)} references.")
        return references


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


class GenericPreprocessor(Preprocessor):

    def __init__(self, config):
        super().__init__(config)
        self._references = self._get_references()


    def __call__(self):
        manager = mp.Manager()
        task_queue = manager.Queue()
        lock = manager.Lock()
        current_transcript_id = manager.Value('i', 1)
        current_read_offset = manager.Value('i', 0)
        processed_transcripts = manager.Value('i', 0)

        logger.info(f"Starting {self._worker_processes} worker processes.")

        workers = [mp.Process(target=self._resquiggle,
                              args=(i,
                                    task_queue,
                                    lock,
                                    current_transcript_id,
                                    current_read_offset,
                                    processed_transcripts,
                                    len(self._references)))
                   for i in range(self._worker_processes)]

        for worker in workers:
            worker.start()

        for ref_id in self._references:
            task_queue.put(ref_id)

        # add poison pills to kill the workers
        for _ in range(self._worker_processes):
            task_queue.put(None)

        for worker in workers:
            worker.join()

        logger.info("Merging tmp databases.")

        # Merge tmp databases
        tmp_dbs = [self._config.get_preprocessing_db() + f".{i}"
                   for i in range(self._worker_processes)]
        self._db.merge_in_databases(tmp_dbs, merge_reads=True)

        logger.info("Creating database indices.")
        self._db.create_indices()


    def _get_references(self):
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "resquiggler preprocessor implementations.")


    def _get_kmer_generator(self, ref_id, ref_seq):
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "resquiggler preprocessor implementations.")


    def _resquiggle(self,
                    idx,
                    task_queue,
                    lock,
                    current_transcript_id,
                    current_read_offset,
                    processed_transcripts,
                    num_transcripts):
        fasta_fh = Fasta(self._config.get_fasta_ref())

        tmp_db_out = PreprocessingDB(self._config.get_preprocessing_db() + f".{idx}",
                                     self._config)
        tmp_db_out.setup()

        with closing(tmp_db_out.connect()) as conn:
            while True:
                ref_id = task_queue.get()
                if ref_id == None:
                    break

                ref_seq = str(fasta_fh[ref_id])
                kmers = list(self._get_kmer_generator(ref_id, ref_seq))
                read_invalid_ratios = self._get_reads_invalid_kmer_ratio(kmers)

                lock.acquire()
                transcript_id = current_transcript_id.value
                current_transcript_id.value += 1
                read_offset = current_read_offset.value
                current_read_offset.value += len(read_invalid_ratios)
                lock.release()

                offsetted_reads = [(read_id, i + read_offset, invalid_ratio)
                                   for i, (read_id, invalid_ratio) in enumerate(read_invalid_ratios.items())]
                kmer_rows = self._process_rows_for_writing(transcript_id, kmers, offsetted_reads)

                with closing(conn.cursor()) as cursor:
                    cursor.execute("begin")
                    cursor.execute(INSERT_TRANSCRIPTS_QUERY, (transcript_id, ref_id))
                    cursor.executemany(INSERT_READS_QUERY, offsetted_reads)
                    cursor.execute("commit")
                PreprocessingDB.write_kmer_rows(conn, kmer_rows)

                lock.acquire()
                processed_transcripts.value += 1
                logger.info(f"Data for reference {processed_transcripts.value}/{num_transcripts} {ref_id} has been saved to the tmp database.")
                lock.release()


    def _process_rows_for_writing(self, transcript_id, kmers_data, offsetted_reads):
        read_to_id = {r[0]: r[1] for r in offsetted_reads}

        processed_kmers = []
        for kmer_data in kmers_data:
            sample_ids = [self._sample_ids[label]
                          for label in kmer_data.sample_labels]
            sample_ids = np.array(sample_ids, dtype=SAMPLE_ID_TYPE)

            read_ids = np.array([read_to_id[read] for read in kmer_data.reads], dtype=READ_ID_TYPE)

            proc_kmer = (transcript_id,
                         kmer_data.pos,
                         kmer_data.kmer,
                         sample_ids.tobytes(),
                         read_ids.tobytes(),
                         kmer_data.intensity.tobytes(),
                         kmer_data.sd.tobytes(),
                         kmer_data.dwell.tobytes())
            processed_kmers.append(proc_kmer)
        return processed_kmers


class Uncalled4Preprocessor(GenericPreprocessor):
    """
    Preprocess a bam file produced by the
    Uncalled4 resquiggler's align command.
    The signal data would be transfered
    to an SQLite DB to be used by Nanocompore
    for later analysis.
    """

    def _get_references(self):
        return self._get_references_from_bams()


    def _get_kmer_generator(self, ref_id, ref_seq):
        uncalled4 = Uncalled4(self._config, ref_id, ref_seq)
        return uncalled4.kmer_data_generator()


class RemoraPreprocessor(GenericPreprocessor):
    """
    Preprocess data using the Remora resquiggler.
    Remora resquiggling is done during the
    preprocessing (as oposed to Uncalled4 and
    Eventalign where the preprocessor will only
    populate the DB with signal alignment data
    that have already been produced by the resquiggler).
    The signal data would be transfered
    to an SQLite DB to be used by Nanocompore
    for later analysis.
    """

    def _get_references(self):
        return self._get_references_from_bams()


    def _get_kmer_generator(self, ref_id, ref_seq):
        remora = Remora(self._config,
                        ref_id=ref_id,
                        start=0,
                        end=len(ref_seq),
                        seq=ref_seq,
                        strand='+')
        return remora.kmer_data_generator()


class EventalignPreprocessor(Preprocessor):
    """
    Takes the output of Nanopolish or F5C eventalign
    command and prepares the data for nanocompore
    by collapsing it and storing it in an SQLite DB
    for later analysis.
    """

    def __init__(self, config):
        super().__init__(config)


    def __call__(self):
        self._reuse_collapsed_files()

        # If any of the samples has a raw eventalign
        # file as input, collapse it.
        tmp_eventalign_dbs = None
        if any('eventalign_db' not in sample_def
               for condition_def in self._config.get_data().values()
               for sample_def in condition_def.values()):
            tmp_eventalign_dbs = self._collapse_eventaligns()

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
            with closing(sqlite3.connect(db_in_path)) as conn_in,\
                 closing(conn_in.cursor()) as db_in,\
                 closing(self._db.connect()) as conn_out,\
                 closing(conn_out.cursor()) as db_out:

                for row in db_in.execute("SELECT name FROM transcripts").fetchall():
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
        tmp_dbs = [self._config.get_preprocessing_db() + f".{i}"
                   for i in range(self._worker_processes)]
        self._db.merge_in_databases(tmp_dbs, merge_reads=False)

        logger.info("Creating database indices.")
        self._db.create_indices()

        if tmp_eventalign_dbs:
            logger.info("Deleting temporary eventalign collapsed databases.")
            for db in tmp_eventalign_dbs:
                Path(db).unlink()


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

        tmp_db_out = PreprocessingDB(self._config.get_preprocessing_db() + f".{idx}",
                                     self._config)
        tmp_db_out.setup()

        while True:
            ref_id = task_queue.get()
            if ref_id == None:
                break

            pos_data = defaultdict(list)
            for sample, db in sample_dbs.items():
                with closing(sqlite3.connect(db)) as conn,\
                     closing(conn.cursor()) as cursor:
                    transcript = cursor.execute("SELECT id FROM transcripts WHERE name = ?",
                                                (ref_id,)).fetchone()
                    if transcript:
                        transcript_id = transcript[0]
                    else:
                        logger.warning(f"Reference {ref_id} not found in the transcripts table for sample {sample}.")
                        continue
                    rows = cursor.execute("SELECT * FROM kmer_data WHERE transcript_id = ?",
                                          (transcript_id,)).fetchall()
                    for row in rows:
                        pos = row[1]
                        pos_data[pos].append((sample, row))

            kmers = self._merge_kmer_data_from_samples(pos_data, sample_read_offsets)
            if not kmers:
                continue

            lock.acquire()
            transcript_id = current_transcript_id.value
            current_transcript_id.value += 1
            lock.release()

            rows = self._process_rows_for_writing(transcript_id, kmers, None)

            with closing(tmp_db_out.connect()) as conn,\
                 closing(conn.cursor()) as cursor:
                cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                               (transcript_id, ref_id))
                PreprocessingDB.write_kmer_rows(conn, rows)
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

            kmers.append(KmerData(None, pos, kmer, samples, read_ids, intensity, sd, dwell, None, self._config))
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

        dbs = []
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
                condition = self._config.sample_to_condition()[sample]
                self._config.get_data()[condition][sample]['eventalign_db'] = db
                dbs.append(db)
        return dbs


    def _get_intermediary_db_name(self, sample):
        kmer_db = self._config.get_preprocessing_db()
        path = Path(kmer_db)
        return str(path.with_name(path.stem + '_' + sample + path.suffix))


    # Since the kmer data we read from the eventalign collapsed db
    # already has internal integer indices instead of uuids,
    # we need a different processing that will not try to map them again.
    # Also, in this case the read mappings are already transfered from
    # the eventalign collapsed db and we don't need to write them here.
    def _process_rows_for_writing(self, transcript_id, kmers_data, read_invalid_ratios=None):
        processed_kmers = []
        for kmer_data in kmers_data:
            sample_ids = [self._sample_ids[label]
                          for label in kmer_data.sample_labels]
            sample_ids = np.array(sample_ids, dtype=SAMPLE_ID_TYPE)

            proc_kmer = (transcript_id,
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
    EventalignCollapser(eventalign, fasta_ref, output, kit, 2)()
    return sample, output

