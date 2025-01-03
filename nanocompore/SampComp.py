import collections, os, traceback, datetime
import sqlite3
import multiprocessing as mp
import time

from contextlib import closing
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor

import torch
import numpy as np

from loguru import logger
from pyfaidx import Fasta

import nanocompore.SampCompResultsmanager as SampCompResultsmanager

from nanocompore.Transcript import Transcript
from nanocompore.TxComp import TxComp
from nanocompore.batch_comp import BatchComp
from nanocompore.common import *
from nanocompore.kmer import KmerData


# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'


DONE_MSG = "done"


class SampComp(object):
    """
    SampComp reads the resquiggled and preprocessed data
    and performs comparison between the samples of the two
    conditions for each position on every whitelisted
    transcript to detect the presence of likely modifications.
    """

    def __init__(self, config, random_state=42):
        """
        Initialise a `SampComp` object and generates a white list of references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        * config
            Config object containing all the parameters for the analysis.
        """
        logger.info("Checking and initialising SampComp")

        # Save init options in dict for later
        log_init_state(loc=locals())

        self._config = config
        self._processing_threads = self._config.get_nthreads() - 1
        self._value_type = get_measurement_type(self._config.get_resquiggler())
        self._sample_labels = self._config.get_sample_labels()
        self._random_state = random_state


    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Starting data processing")

        resultsManager = SampCompResultsmanager.resultsManager(self._config)

        logger.info(f"Starting {self._processing_threads} worker processes.")

        manager = mp.Manager()
        task_queue = manager.Queue()
        result_queue = manager.Queue()
        if self._config.get_device() == 'cpu':
            logger.info("Will process transcripts using the CPU.")
            transcript_processor = self._process_transcripts
        else:
            logger.info("Will process transcripts using the GPU.")
            transcript_processor = self._process_transcripts_gpu
        workers = [mp.Process(target=transcript_processor,
                              args=(i, task_queue, result_queue))
                   for i in range(self._processing_threads)]

        for worker in workers:
            worker.start()

        transcripts = self._get_transcripts_for_processing()
        if self._config.get_result_exists_strategy() == 'continue':
            all_transcripts_num = len(transcripts)
            transcripts = resultsManager.filter_already_processed_transcripts(transcripts)
            already_processed = all_transcripts_num - len(transcripts)
            logger.info(f"Found a total of {all_transcripts_num}. {already_processed} have already been processed in previous runs. Will process only for the remaining {len(transcripts)}.")
        else:
            logger.info(f"Found a total of {len(transcripts)} transcripts for processing.")

        
        # for ref_ids in chunks(transcripts, 1):
        #     task_queue.put(ref_ids)
        for transcript_ref in transcripts:
            # if "ENST00000274606.8" in transcript_ref.ref_id:
            task_queue.put(transcript_ref)

        # add poison pills to kill the workers
        for _ in range(self._processing_threads):
            task_queue.put(None)

        finished_workers = 0
        while True:
            msg = result_queue.get()
            if msg == DONE_MSG:
                finished_workers += 1
                if finished_workers == self._processing_threads:
                    break
                continue

            transcript, results = msg
            logger.info(f"Main process received results for {transcript.name}. Will save them now.")
            resultsManager.save_results(transcript, results)

        for worker in workers:
            worker.join()

        resultsManager.finish()


    def _process_transcripts(self, worker_id, task_queue, result_queue):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        batch_comp = BatchComp(config=self._config)
        kit = self._config.get_kit()

        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            while True:
                transcript_ref = task_queue.get()
                if transcript_ref is None:
                    logger.info(f"Worker {worker_id} has finished and will terminate.")
                    result_queue.put(DONE_MSG)
                    return

                # tx_id_to_transcript = {ref.id: Transcript(ref_id=ref.ref_id,
                #                                           ref_seq=str(fasta_fh[ref.ref_id]))
                #                        for ref in transcript_refs}
                transcript = Transcript(id=transcript_ref.id,
                                        ref_id=transcript_ref.ref_id,
                                        ref_seq=str(fasta_fh[transcript_ref.ref_id]))

                try:
                    kmer_data_list = self._read_transcript_kmer_data([transcript_ref], cursor)
                    transcript, results = batch_comp.compare_transcripts(transcript, kmer_data_list)

                    seq = transcript.seq
                    kmers = results.pos.apply(lambda pos: get_pos_kmer(pos, seq, kit))
                    results['kmer'] = kmers.apply(encode_kmer)
                    result_queue.put((transcript, results))
                    # for tx_id in all_test_results['transcript_id'].unique():
                    #     df = all_test_results[all_test_results.transcript_id == tx_id].copy()
                    #     transcript = tx_id_to_transcript[tx_id]
                    #     ref_id = transcript.name
                    #     seq = transcript.seq
                    #     kmers = df.pos.apply(lambda pos: get_pos_kmer(pos, seq, kit))
                    #     df['kmer'] = kmers.apply(encode_kmer)
                    #     result_queue.put((ref_id, df))

                except Exception:
                    msg = traceback.format_exc()
                    logger.error(f"Error in Worker {worker_id}: {msg}")


    def _process_transcripts_gpu(self, worker_id, task_queue, result_queue):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        batch_comp = BatchComp(config=self._config)
        kit = self._config.get_kit()

        manager = mp.Manager()
        max_threads = manager.Value('i', 1)
        threads = mp.pool.ThreadPool(max_threads.value)
        active_threads = manager.Value('i', 0)
        retry_queue = manager.Queue()

        # TODO add lock
        async_results = {}
        mem_usage_history = []

        # kill_proc = False
        input_finished = False

        def callback(args):
            transcript, results = args
            logger.info(f"Transcript {transcript.name} has been processed, will send the results for saving.")
            try:
                active_threads.value -= 1
                del async_results[transcript.name]
                # ref_id = transcript.name
                seq = transcript.seq
                kmers = results.pos.apply(lambda pos: get_pos_kmer(pos, seq, kit))
                results['kmer'] = kmers.apply(encode_kmer)
                result_queue.put((transcript, results))
                logger.info(f"Transcript {transcript.name} has been processed and was sent for saving.")
            except Exception as e:
                msg = traceback.format_exc()
                logger.error(f"Got an exception in the async callback: {msg}")

        def error_callback(exception):
            transcript = exception.transcript
            logger.error(f"Worker {worker_id} got exception {exception} while processing {transcript.name}")
            active_threads.value -= 1
            del async_results[transcript.name]
            if isinstance(exception, torch.OutOfMemoryError):
                logger.error(f"Worker {worker_id} got OutOfMemory error while processing {transcript.name}. Will add it to the retry queue.")
                # if max_threads.value > 1:
                #     max_threads.value -= 1
                max_threads.value = max(active_threads.value - 1, 1)
                retry_queue.put(TranscriptRow(transcript.name, transcript.id))
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            else:
                # TODO: format stack trace and log that
                # msg = traceback.format_exception(exception)
                logger.error(f"Error in Worker {worker_id}: {exception}")

        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            t = time.time()
            while True:
                free, total = torch.cuda.mem_get_info()
                mem_usage = (total - free) / total
                if time.time() - t > 10:
                    mem_usage_history.insert(0, mem_usage)
                    if len(mem_usage_history) > 3:
                        mem_usage_history.pop()
                    if len(mem_usage_history) == 3 and np.mean(mem_usage_history) < 0.65:
                        # max_threads.value += 1
                        mem_usage_history = []
                    logger.info(f"Worker {worker_id}: Active threads: {active_threads.value}, max threads: {max_threads.value}, mem_usage: {mem_usage}")
                    t = time.time()

                # if kill_proc and all(async_results.ready() for async_res in async_results.values()):
                if input_finished and len(async_results) == 0 and retry_queue.empty():
                    logger.info(f"Worker {worker_id} has finished and will terminate.")
                    threads.close()
                    result_queue.put(DONE_MSG)
                    return
                
                if active_threads.value < max_threads.value and mem_usage < 0.80:
                    if not input_finished:
                        transcript_ref = task_queue.get()
                        if transcript_ref is None:
                            # threads.close()
                            # kill_proc = True
                            input_finished = True
                            continue
                    elif not retry_queue.empty():
                        transcript_ref = retry_queue.get()
                    else:
                        continue


                    transcript = Transcript(id=transcript_ref.id,
                                            ref_id=transcript_ref.ref_id,
                                            ref_seq=str(fasta_fh[transcript_ref.ref_id]))
                    kmer_data_list = self._read_transcript_kmer_data([transcript_ref], cursor)
                    logger.info(f"Worker {worker_id} submits transcript {transcript_ref.ref_id} for processing to the thread pool.")
                    async_res = threads.apply_async(batch_comp.compare_transcripts,
                                                    args=(transcript, kmer_data_list,),
                                                    callback=callback,
                                                    error_callback=error_callback)
                    async_results[transcript.name] = async_res
                    active_threads.value += 1
                else:
                    time.sleep(1)


                # transcript_refs = task_queue.get()
                # if transcript_refs is None:
                #     logger.info(f"Worker {worker_id} has finished and will terminate.")
                #     result_queue.put(DONE_MSG)
                #     return

                # tx_id_to_transcript = {ref.id: Transcript(ref_id=ref.ref_id,
                #                                           ref_seq=str(fasta_fh[ref.ref_id]))
                #                        for ref in transcript_refs}

                # try:
                #     kmer_data_list = self._read_transcript_kmer_data(transcript_refs, cursor)
                #     all_test_results = batch_comp.compare_transcripts(kmer_data_list)

                #     for tx_id in all_test_results['transcript_id'].unique():
                #         df = all_test_results[all_test_results.transcript_id == tx_id].copy()
                #         transcript = tx_id_to_transcript[tx_id]
                #         ref_id = transcript.name
                #         seq = transcript.seq
                #         kmers = df.pos.apply(lambda pos: get_pos_kmer(pos, seq, kit))
                #         df['kmer'] = kmers.apply(encode_kmer)
                #         result_queue.put((ref_id, df))

                # except Exception:
                #     msg = traceback.format_exc()
                #     logger.error(f"Error in Worker {worker_id}: {msg}")


    def _read_transcript_kmer_data(self, transcript_refs, cursor):
        params = ','.join(['?' for _ in transcript_refs])
        res = cursor.execute(f"""SELECT *
                                 FROM kmer_data
                                 WHERE transcript_id IN ({params})""",
                             [ref.id for ref in transcript_refs])

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
        # discard excesive reads
        max_reads = self._config.get_downsample_high_coverage()
        return [kmer.subsample_reads(max_reads, random_state=self._random_state)
                for kmer in kmers]


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
        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            query = f"""
            SELECT DISTINCT t.reference, t.id
            FROM kmer_data kd
            INNER JOIN transcripts t ON kd.transcript_id = t.id
            """
            return {TranscriptRow(row[0], row[1])
                    for row in cursor.execute(query).fetchall()}

