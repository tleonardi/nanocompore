# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
import collections, os, traceback, datetime
import multiprocessing as mp

# Third party
from loguru import logger
import yaml
#from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta

# Local package
from nanocompore.common import *
import nanocompore.SampCompResultsmanager as SampCompResultsmanager
import nanocompore.Transcript as Transcript
import nanocompore.TxComp as TxComp
import nanocompore.Whitelist as Whitelist
import nanocompore.Experiment as Experiment

import nanocompore as pkg

from concurrent.futures import ProcessPoolExecutor, as_completed


# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampComp(object):
    """ Init analysis and check args"""

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

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

        # Set private args
        self._experiment = Experiment.Experiment(config)

        logger.info("Starting to whitelist the reference IDs")
        self._Whitelist = Whitelist.Whitelist(self._experiment, self._config)

        self._valid_transcripts = self._Whitelist.ref_id_list
        self._processing_threads = self._config.get_nthreads() - 1


    #~~~~~~~~~~~~~~Call METHOD~~~~~~~~~~~~~~#
    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Starting data processing")

        resultsManager = SampCompResultsmanager.resultsManager(self._config)
        try:
            with ProcessPoolExecutor(max_workers=self._processing_threads) as executor:
                futures = [executor.submit(self._processTx, ref_id)
                           for ref_id in self._valid_transcripts]
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


    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#


    def _processTx(self, ref_id):
        logger.debug(f"Worker thread processing new item from in_q: {ref_id}")
        fasta_fh = Fasta(self._config.get_fasta_ref())
        ref_seq = str(fasta_fh[ref_id])
        txComp = TxComp.TxComp(experiment=self._experiment,
                               config=self._config,
                               random_state = 26)
        transcript = Transcript.Transcript(ref_id=ref_id,
                                           experiment=self._experiment,
                                           ref_seq=ref_seq,
                                           config=self._config)
        try:
            results = txComp.txCompare(transcript)
            return ref_id, results
        except Exception as e:
            import traceback
            traceback.print_stack(e)
            logger.error(f"Error in Worker for {ref_id}: {e}")


    def _writeResults(self, out_q, error_q):
        #TODO need to determine if this is necessary
        #self.resultsManager.saveExperimentMetadata()

        try:
            n_tx = 0
            n_pos = 0
            for _ in range(self._processing_threads):
                for tx, result in iter(out_q.get, None):
                    if result:
                        logger.debug(f"Writer thread adding results data from {tx}")
                        n_tx += 1
                        n_pos = len([x for x in result if type(x) == int])
                        self.resultsManager.saveData(tx, result, self._config)
            self.resultsManager.finish()
        except:
            logger.error("Error writing results to database")
            error_q.put(NanocomporeError(traceback.format_exc()))
        finally:
            logger.debug("Written Transcripts:{} Valid positions:{}".format(n_tx, n_pos))
            self.resultsManager.closeDB()
            logger.info ("All Done. Transcripts processed: {}".format(n_tx))
            # Kill error queue with poison pill
            error_q.put(None)

