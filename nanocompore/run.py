import multiprocessing as mp
import os
import queue
import sqlite3
import time
import traceback

from contextlib import closing

import numpy as np
import torch

from loguru import logger
from pyfaidx import Fasta

from nanocompore.common import Counter
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.common import TranscriptRow
from nanocompore.common import encode_kmer
from nanocompore.common import get_measurement_type
from nanocompore.common import get_pos_kmer
from nanocompore.common import log_init_state
from nanocompore.comparisons import TranscriptComparator
from nanocompore.database import ResultsDB
from nanocompore.kmer import KmerData
from nanocompore.postprocessing import Postprocessor
from nanocompore.transcript import Transcript


class RunCmd(object):
    """
    RunCmd reads the resquiggled and preprocessed data
    and performs comparison between the samples of the two
    conditions for each position on every transcript to
    detect the presence of likely modifications.
    """

    def __init__(self, config, random_state=42):
        """
        Initialise a `RunCmd` object and generates a white list of references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        * config
            Config object containing all the parameters for the analysis.
        """
        logger.info("Checking and initialising RunCmd")

        # Save init options in dict for later
        log_init_state(loc=locals())

        self._config = config
        self._workers_with_device = self._get_workers_with_device()
        self._value_type = get_measurement_type(self._config.get_resquiggler())
        self._sample_labels = self._config.get_sample_labels()
        self._random_state = random_state


    def __call__(self):
        """
        The main method that will perform the comparative analysis.
        The function will get transcripts for processing from the input
        database, start worker processes to compare them and write the results
        in the output database, and finally will perform the postprocessing
        to export TSVs with the output.
        """
        logger.info("Starting data processing")

        if self._config.get_progress():
            print("Getting the list of transcripts for processing...")
        transcripts = self._get_transcripts_for_processing()
        db_path = ResultsDB(self._config, init_db=False).db_path
        db_exists = os.path.isfile(db_path)
        if db_exists and self._config.get_result_exists_strategy() == 'continue':

            all_transcripts_num = len(transcripts)
            transcripts = self._filter_processed_transcripts(transcripts)
            already_processed = all_transcripts_num - len(transcripts)
            logger.info(f"Found a total of {all_transcripts_num}. " + \
                        f"{already_processed} have already been processed in previous runs. " + \
                        f"Will process only for the remaining {len(transcripts)}.")
        else:
            # Initialize the database
            ResultsDB(self._config, init_db=True)
            logger.info(f"Found a total of {len(transcripts)} transcripts for processing.")

        task_queue = mp.JoinableQueue()
        for transcript_ref in transcripts:
            task_queue.put((transcript_ref, 0))
        logger.info("All transcripts have been sent to the queue")

        if self._config.get_progress():
            print("Starting to compare the conditions for all transcripts...")
        logger.info(f"Starting {len(self._workers_with_device)} worker processes.")
        db_lock = mp.Lock()
        sync_lock = mp.Lock()
        manager = mp.Manager()
        num_finished = manager.Value('i', 0)
        workers = [mp.Process(target=self._worker_process,
                              args=(worker, task_queue, db_lock, sync_lock, num_finished, device))
                   for worker, device in self._workers_with_device.items()]

        for worker in workers:
            worker.start()

        # helper variables for the progress bar
        max_len = 0
        start_time = time.time()

        # Wait for all workers to finish while printing
        # the progress bar if it was requested.
        # We use this slightly weird pattern of waiting
        # for processes instead of a standard loop with
        # unbounded join, because otherwise we'll have
        # to start a separate thread for the progress bar
        # monitoring and that doesn't work well with the
        # multiprocess Manager.
        while len(workers) > 0:
            for i, worker in enumerate(workers):
                # Try to join the worker with 1s timeout
                worker.join(1)
                # If the worker has completed we remove it from
                # the worker pool.
                if worker.exitcode is not None:
                    workers.pop(i)
                if self._config.get_progress():
                    with sync_lock:
                        num_done = num_finished.value
                    self._print_progress_bar(num_done, len(transcripts), max_len, start_time)

        if self._config.get_progress():
            print("\nAll transcripts have been processed.")
            print("Starting postprocessing...")

        postprocessor = Postprocessor(self._config)
        postprocessor()

        if self._config.get_progress():
            print("Done.")


    def _worker_process(self, worker_id, task_queue, db_lock, sync_lock, num_finished, device):
        """
        The main worker process function. It will
        continuously take tasks (transcripts) from the
        task queue, process them and write the results
        to the database. Will quit when the task queue
        is empty.

        :param worker_id int: Id of the worker
        :param task_queue multiprocesing.JoinableQueue: Queue with tasks for processing.
        :param db_lock multiprocessing.Lock: Lock to prevent parallel writing to the SQLite DB.
        :param sync_lock multiprocessing.Lock: Lock for synchronisation.
        :param num_finished multiprocessing.managers.SyncManager.Value: number of processed transcripts.
        :param device str: which device to use for computation (e.g. cpu or gpu).
        """
        db_manager = ResultsDB(self._config, init_db=False)
        fasta_fh = Fasta(self._config.get_fasta_ref())
        comparator = TranscriptComparator(config=self._config)
        kit = self._config.get_kit()
        logger.info(f"Worker {worker_id} started with pid {os.getpid()}.")

        db = self._config.get_preprocessing_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            while True:
                # It's weird that we have to use a lock here,
                # but for some reason calling get in parallel from
                # multiple processes can cause a queue.Empty
                # exception to be raised even for non-empty queues...
                # https://github.com/python/cpython/issues/87302
                with sync_lock:
                    try:
                        logger.info(f"Worker {worker_id} getting a task from a queue with size {task_queue.qsize()}")
                        msg = task_queue.get(block=False)
                    except queue.Empty:
                        logger.info(f"Worker {worker_id} could not find more tasks in the queue and will terminate.")
                        break

                try:
                    transcript_ref, retries = msg
                    if retries > 2:
                        logger.error(f"Worker {worker_id}: Transcript {transcript_ref.ref_id} was retried too many times. Won't retry again.")
                        continue
                    logger.info(f"Worker {worker_id} starts processing {transcript_ref.ref_id}")
                    transcript = Transcript(id=transcript_ref.id,
                                            ref_id=transcript_ref.ref_id,
                                            ref_seq=str(fasta_fh[transcript_ref.ref_id]))

                    kmer_data_list = self._read_transcript_kmer_data(transcript_ref, cursor)
                    transcript, results = comparator.compare_transcript(transcript,
                                                                        kmer_data_list,
                                                                        device)
                    if results is None:
                        continue

                    seq = transcript.seq
                    kmers = results.pos.apply(lambda pos: get_pos_kmer(pos, seq, kit))
                    results['kmer'] = kmers.apply(encode_kmer)
                    with db_lock:
                        logger.info(f"Worker {worker_id}: Saving the results for {transcript.name}")
                        db_manager.save_test_results(transcript, results)
                except torch.OutOfMemoryError:
                    task = TranscriptRow(transcript.name, transcript.id)
                    logger.info(f"Worker {worker_id} returns {transcript_ref.ref_id} to " + \
                                 "the queue for later processing because it failed with OutOfMemory.")
                    with sync_lock:
                        task_queue.put((task, retries + 1))
                except:
                    msg = traceback.format_exc()
                    logger.error(f"Error in Worker {worker_id}: {msg}")
                finally:
                    logger.info(f"Worker {worker_id} is finishing a task.")
                    with sync_lock:
                        num_finished.value += 1
                        task_queue.task_done()


    def _read_transcript_kmer_data(self, transcript_ref, cursor):
        res = cursor.execute("""SELECT *
                                FROM kmer_data
                                WHERE transcript_id = ?""",
                             (transcript_ref.id,))

        kmers = [self._db_row_to_kmer(row)
                 for row in res.fetchall()]

        # After reading the data for all kmers,
        # we can collect all reads and select
        # only the ones that pass the filtering
        # criteria.
        all_read_ids = {read_id
                        for kmer in kmers
                        for read_id in kmer.reads}
        valid_read_ids = self._get_valid_read_ids(cursor, all_read_ids)
        logger.debug(f"Found {len(valid_read_ids)} valid reads out of {len(all_read_ids)}.")
        # discard invalid reads
        kmers = [self._filter_reads(kmer, valid_read_ids) for kmer in kmers]
        # get only kmers that have the minimum number of reads required
        kmers = [kmer for kmer in kmers if self._enough_reads_in_kmer(kmer)]
        return kmers


    def _filter_reads(self, kmer, valid_read_ids):
        valid_mask = [read in valid_read_ids for read in kmer.reads]
        return KmerData(kmer.transcript_id,
                        kmer.pos,
                        kmer.kmer,
                        kmer.sample_labels[valid_mask],
                        kmer.reads[valid_mask],
                        kmer.intensity[valid_mask],
                        kmer.sd[valid_mask],
                        kmer.dwell[valid_mask],
                        None, # We don't have validity data here
                        self._config)


    def _enough_reads_in_kmer(self, kmer):
        condition_counts = Counter(kmer.condition_labels)
        return all(condition_counts.get(cond, 0) >= self._config.get_min_coverage()
                   for cond in self._config.get_condition_labels())


    def _db_row_to_kmer(self, row):
        read_ids = np.frombuffer(row[4], dtype=READ_ID_TYPE)
        samples = np.array([self._sample_labels[i]
                            for i in np.frombuffer(row[3], dtype=SAMPLE_ID_TYPE)])
        intensity = np.frombuffer(row[5], dtype=self._value_type)
        sd = np.frombuffer(row[6], dtype=self._value_type)
        dwell = np.frombuffer(row[7], dtype=self._value_type)

        return KmerData(row[0], # transcript id
                        row[1], # pos
                        row[2], # kmer
                        samples,
                        read_ids,
                        intensity,
                        sd,
                        dwell,
                        None, # we don't need validity data here
                        self._config)


    def _get_valid_read_ids(self, cursor, read_ids):
        get_ids_query = """
        SELECT id
        FROM reads
        WHERE id IN ({})
          AND invalid_kmers <= ?
        """.format(','.join(map(str, read_ids)))
        res = cursor.execute(get_ids_query,
                             (self._config.get_max_invalid_kmers_freq(),))

        return {row[0] for row in res.fetchall()}


    def _get_transcripts_for_processing(self):
        db = self._config.get_preprocessing_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            query = """
            SELECT DISTINCT t.name, t.id
            FROM kmer_data kd
            INNER JOIN transcripts t ON kd.transcript_id = t.id
            """
            return {TranscriptRow(row[0], row[1])
                    for row in cursor.execute(query).fetchall()}


    def _filter_processed_transcripts(self, transcripts):
        db = ResultsDB(self._config, init_db=False)
        logger.info(f'Preprocessing db: {db._db_path}')
        existing_transcripts = set(db.get_transcripts())
        return [tx for tx in transcripts if tx.ref_id not in existing_transcripts]


    def _get_workers_with_device(self):
        devices = self._config.get_devices()
        result = {}
        if isinstance(devices, dict):
            worker = 0
            for device, procs in devices.items():
                for _ in range(procs):
                    result[worker] = device
                    worker += 1
        else:
            if not isinstance(devices, list):
                devices = [devices]
            for worker in range(self._config.get_nthreads() - 1):
                result[worker] = devices[worker % len(devices)]
        return result


    def _print_progress_bar(self, num_done, total_transcripts, max_len, start_time):
        progress_bar_len = 30
        perc_done = num_done / total_transcripts
        done_chars = int(perc_done*progress_bar_len)
        done_bar = '▮' * done_chars
        remaining_bar = '-' * (progress_bar_len - done_chars)
        running_time = time.time() - start_time
        if num_done == 0:
            eta = 'ETA: N/A'
        else:
            total_expected = running_time * total_transcripts / num_done
            remaining_sec = total_expected - running_time
            hours, minutes, _ = self._format_time(remaining_sec)
            eta = f'ETA: {hours:02d}h:{minutes:02d}m remaining'
        hours, minutes, seconds = self._format_time(running_time)
        running = f'⏱  {hours:02d}h:{minutes:02d}m:{seconds:02d}s'
        speed = num_done / (running_time / 60)
        bar = f'{running} |{done_bar}{remaining_bar}| ' + \
              f'{num_done:,}/{total_transcripts:,} ' + \
              f'({perc_done*100:3.1f}%) [{speed:.1f}tx/m] | {eta}'
        if len(bar) > max_len:
            max_len = len(bar)
        if len(bar) < max_len:
            bar = f'{bar}{(max_len - len(bar)) * " "}'
        print(bar, end='\r', flush=True)
        return max_len


    def _format_time(self, seconds):
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return int(hours), int(minutes), int(seconds)

