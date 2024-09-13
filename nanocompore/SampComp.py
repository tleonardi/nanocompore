import collections, os, traceback, datetime
import sqlite3
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

        try:
            with ProcessPoolExecutor(max_workers=self._processing_threads) as executor:
                futures = [executor.submit(self._processTx, tx)
                           for tx in self._get_transcripts_for_processing()]
                for future in as_completed(futures):
                    ref_id, results = future.result()
                    if len(result) == 0:
                        continue
                    try:
                        resultsManager.saveData(ref_id, results, self._config)
                    except Exception:
                        logger.error("Error writing results to database")
                        raise NanocomporeError(traceback.format_exc())
            resultsManager.finish()
        finally:
            resultsManager.closeDB()


    def _processTx(self, transcript_ref):
        logger.debug(f"Worker thread starts processing new transcript: {transcript_ref.ref_id}")
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[transcript_ref.ref_id])
        txComp = TxComp(experiment=self._experiment,
                        config=self._config,
                        random_state=26)
        transcript = Transcript(ref_id=transcript_ref.ref_id,
                                experiment=self._experiment,
                                ref_seq=ref_seq,
                                config=self._config)

        kmer_data_list = self._read_transcript_kmer_data(transcript_ref)
        try:
            results = txComp.txCompare(kmer_data_list, transcript)
            return transcript_ref.ref_id, results
        except Exception as e:
            traceback.print_stack(e)
            logger.error(f"Error in Worker for {ref_id}: {e}")


    def _read_transcript_kmer_data(self, transcript_ref):
        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
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
            return [TranscriptRow(row[0], row[1])
                    for row in cursor.execute("SELECT reference, id FROM transcripts").fetchall()]

