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
        self._processing_threads = self._config.get_nthreads() - 2


    #~~~~~~~~~~~~~~Call METHOD~~~~~~~~~~~~~~#
    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Starting data processing")

        in_q = mp.Queue(maxsize = 100)
        out_q = mp.Queue(maxsize = 100)
        error_q = mp.Queue(maxsize = 100)

        for tx in self._valid_transcripts:
            in_q.put(tx)
        
        self.resultsManager = SampCompResultsmanager.resultsManager(self._config)

        self.txComp = TxComp.TxComp(experiment=self._experiment,
                                    config=self._config,
                                    random_state = 26)
        # Define processes
        processes = list()
        processes.append(mp.Process(target=self._writeResults, args=(out_q, error_q)))
        for _ in range(self._processing_threads):
            in_q.put(None)
            processes.append(mp.Process(target=self._processTx, args=(in_q, out_q, error_q)))

        try:
            #Start all processes
            self._fasta_fh = Fasta(self._config.get_fasta_ref())
            for p in processes:
                p.start()

            # Monitor error queue
            for tb in iter(error_q.get, None):
                logger.trace("Error caught from error_q")
                raise NanocomporeError(tb)
            logger.debug("Error queue was closed")
            
            # Soft processes and queues stopping
            logger.debug("Waiting for all processes to be joined")
            for p in processes:
                p.join()
            logger.debug("All processes joined successfully")

            logger.debug("Closing all queues")
            for q in (in_q, out_q, error_q):
                q.close()
            logger.debug("All queues were closed")

            self.resultsManager.finish()

        except Exception as E:
            logger.error("An error occured. Killing all processes and closing queues\n")
            try:
                for p in processes:
                    p.terminate()
                for q in (in_q, out_q, error_q):
                    q.close()
            except:
                logger.error("An error occured while trying to kill processes\n")
            raise E
        finally:
            self.resultsManager.closeDB()
            self._fasta_fh.close()


    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#
    def _processTx(self, in_q, out_q, error_q):
        logger.debug("Worker thread started")
        try:
            n_tx = 0
            for ref_id in iter(in_q.get, None):
                logger.debug(f"Worker thread processing new item from in_q: {ref_id}")
                ref_seq = str(self._fasta_fh[ref_id])
                transcript = Transcript.Transcript(ref_id=ref_id,
                                                   experiment=self._experiment,
                                                   ref_seq=ref_seq,
                                                   config=self._config)
                try:
                    n_tx += 1
                    logger.debug(f"Collecting data for {ref_id}")
                    results = self.txComp.txCompare(transcript)
                    out_q.put((transcript.name, results))
                except:
                    logger.debug(f"Insufficent coverage for {ref_id} skipping transcript")

        except:
            logger.error("Error in Worker")
            error_q.put(NanocomporeError(traceback.format_exc()))

        # Deal poison pill and close file pointer
        finally:
            logger.debug(f"Processed Transcrits: {n_tx}")
            logger.debug("Adding poison pill to out_q")
            out_q.put(None)

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
                        self.resultsManager.saveData(tx, result)
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

    #TODO delete?
    def _list_refid(self, in_q, error_q):
        """Add valid refid from whitelist to input queue to dispatch the data among the workers"""
        n_tx = 0
        try:
            for tx in self._valid_transcripts:
                logger.debug(f"Adding {tx} to in_q")
                in_q.put((tx))
                n_tx+=1

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in Reader")
            error_q.put(traceback.format_exc())

        # Deal poison pills
        finally:
            for i in range (self._nthreads):
                in_q.put(None)
            logger.debug(f"Parsed transcripts:{n_tx}")