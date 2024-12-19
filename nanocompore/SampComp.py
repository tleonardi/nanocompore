import collections, os, traceback, datetime
import sqlite3
import multiprocessing as mp
from contextlib import closing
from concurrent.futures import ProcessPoolExecutor, as_completed

from loguru import logger
#from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta

from nanocompore.common import *
import nanocompore.SampCompResultsmanager as SampCompResultsmanager
from nanocompore.Transcript import Transcript
from nanocompore.TxComp import TxComp
from nanocompore.Experiment import Experiment
from nanocompore.kmer import KmerData


# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

class SampComp(object):
    """
    SampComp reads the resquiggled and preprocessed data
    and performs comparison between the samples of the two
    conditions for each position on every whitelisted
    transcript to detect the presence of likely modifications.
    """

    def __init__(self, config):
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
        self._experiment = Experiment(config)
        self._processing_threads = self._config.get_nthreads() - 1
        self._value_type = get_measurement_type(self._config.get_resquiggler())
        self._sample_labels = self._experiment.get_sample_labels()


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
        workers = [mp.Process(target=self._processTx,
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
        for ref_id in transcripts:
            task_queue.put(ref_id)

        # add poison pills to kill the workers
        for _ in range(self._processing_threads):
            task_queue.put(None)

        try:
            finished_workers = 0
            while True:
                msg = result_queue.get()
                if msg == "done":
                    finished_workers += 1
                    if finished_workers == self._processing_threads:
                        break
                    continue
                ref_id, results = msg
                logger.info(f"Got {type(results)}")
                logger.info(f"Got {len(results)} results for saving for {ref_id}")
                if len(results) == 0:
                    continue
                logger.info(f"Starts saving results for {ref_id} to the SampComp database.")
                resultsManager.saveData(ref_id, results, self._config)
                logger.info(f"Saved results for {ref_id} to the SampComp database.")

            for worker in workers:
                worker.join()

            resultsManager.finish()
        finally:
            resultsManager.closeDB()


    def _processTx(self, i, task_queue, result_queue):
        fasta_fh = Fasta(self._config.get_fasta_ref())
        txComp = TxComp(experiment=self._experiment,
                        config=self._config,
                        random_state=26)

        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            while True:
                transcript_ref = task_queue.get()
                if transcript_ref is None:
                    logger.debug(f"Worker {i} has finished and will terminate.")
                    result_queue.put("done")
                    break

                try:
                    logger.debug(f"Worker {i} got new transcript for processing: {transcript_ref.ref_id}")
                    ref_seq = str(fasta_fh[transcript_ref.ref_id])
                    transcript = Transcript(ref_id=transcript_ref.ref_id,
                                            experiment=self._experiment,
                                            ref_seq=ref_seq,
                                            config=self._config)
                    logger.debug(f"Worker thread starts reading data for transcript: {transcript_ref.ref_id}")
                    kmer_data_list = self._read_transcript_kmer_data(transcript_ref, cursor)
                    logger.debug(f"Worker thread starts comparing data for transcript: {transcript_ref.ref_id}")
                    results = txComp.txCompare(kmer_data_list, transcript)
                    logger.debug(f"Worker thread finished comparing data for transcript: {transcript_ref.ref_id}. Results: {len(results)}")
                except Exception as e:
                    traceback.print_stack(e)
                    logger.error(f"Error in Worker for {ref_id}: {e}")


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
        logger.debug(f"Found {len(valid_read_ids)} valid reads out of {len(all_read_ids)} total for {transcript_ref.ref_id}.")
        return [self._filter_reads(kmer, valid_read_ids) for kmer in kmers]


    def _filter_reads(self, kmer, valid_read_ids):
        valid_mask = [read in valid_read_ids for read in kmer.reads]
        return KmerData(kmer.pos,
                        kmer.kmer,
                        kmer.sample_labels[valid_mask],
                        None, # read ids are not necessary for sampcomp
                        kmer.intensity[valid_mask],
                        kmer.sd[valid_mask],
                        kmer.dwell[valid_mask],
                        None, # We don't have validity data here
                        self._experiment)


    def _db_row_to_kmer(self, row):
        read_ids = np.frombuffer(row[4], dtype=READ_ID_TYPE)
        samples = np.array([self._sample_labels[i]
                            for i in np.frombuffer(row[3], dtype=SAMPLE_ID_TYPE)])
        intensity = np.frombuffer(row[5], dtype=self._value_type)
        sd = np.frombuffer(row[6], dtype=self._value_type)
        dwell = np.frombuffer(row[7], dtype=self._value_type)

        return KmerData(row[1], # pos
                        row[2], # kmer
                        samples,
                        read_ids,
                        intensity,
                        sd,
                        dwell,
                        None, # we don't need validity data here
                        self._experiment)


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
            SELECT t.reference, t.id
            FROM kmer_data kd
            INNER JOIN transcripts t ON kd.transcript_id = t.id
            WHERE num_reads >= {self._config.get_min_coverage()}
            """
            return {TranscriptRow(row[0], row[1])
                    for row in cursor.execute(query).fetchall()}

