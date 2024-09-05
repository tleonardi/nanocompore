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

import nanocompore as pkg


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


    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Starting data processing")

        resultsManager = SampCompResultsmanager.resultsManager(self._config)

        try:
            with ProcessPoolExecutor(max_workers=self._processing_threads) as executor:
                futures = [executor.submit(self._processTx, ref_id)
                           for ref_id in self._get_transcripts_for_processing()]
                for future in as_completed(futures):
                    ref_id, results = future.result()
                    try:
                        resultsManager.saveData(ref_id, results, self._config)
                    except Exception:
                        logger.error("Error writing results to database")
                        raise NanocomporeError(traceback.format_exc())
            resultsManager.finish()
        finally:
            resultsManager.closeDB()


    def _processTx(self, ref_id):
        logger.debug(f"Worker thread processing new item from in_q: {ref_id}")
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])
        txComp = TxComp(experiment=self._experiment,
                        config=self._config,
                        random_state=26)
        transcript = Transcript(ref_id=ref_id,
                                experiment=self._experiment,
                                ref_seq=ref_seq,
                                config=self._config)

        kmer_data_list = self._read_transcript_kmer_data(ref_id)
        try:
            results = txComp.txCompare(kmer_data_list, transcript)
            return ref_id, results
        except Exception as e:
            traceback.print_stack(e)
            logger.error(f"Error in Worker for {ref_id}: {e}")


    def _read_transcript_kmer_data(self, ref_id):
        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            res = cursor.execute("SELECT * FROM kmer_data WHERE transcript = ?", (ref_id,))
            sample_labels = self._experiment.get_sample_labels()
            value_type = get_measurement_type(self._config.get_resquiggler())
            return [KmerData(row[1], # pos
                             row[2], # kmer
                             [sample_labels[i]
                              for i in np.frombuffer(row[3], dtype=SAMPLE_ID_TYPE)],
                             None, # read ids are not necessary for sampcomp
                             np.frombuffer(row[5], dtype=value_type), # intensity
                             np.frombuffer(row[6], dtype=value_type), # sd
                             np.frombuffer(row[7], dtype=value_type), # dwell
                             self._experiment)
                    for row in res.fetchall()]


    def _get_transcripts_for_processing(self):
        db = self._config.get_kmer_data_db()
        with closing(sqlite3.connect(db)) as conn,\
             closing(conn.cursor()) as cursor:
            return [row[0]
                    for row in cursor.execute("SELECT DISTINCT transcript FROM kmer_data").fetchall()]

