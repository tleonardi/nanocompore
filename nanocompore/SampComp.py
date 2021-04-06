# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import shelve
import multiprocessing as mp
import traceback
import datetime
import os

# Third party
from loguru import logger
from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta

# Local package
from nanocompore.common import *
from nanocompore.DataStore import *
from nanocompore.Whitelist import Whitelist
from nanocompore.TxComp import txCompare
from nanocompore.SampCompDB import SampCompDB
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

    def __init__(self,
                 input_db_path:str,
                 output_db_path:str,
                 sample_dict:dict,
                 fasta_fn:str = "",
                 overwrite:bool = False,
                 whitelist:Whitelist = None,
                 comparison_methods:list = ["GMM", "KS"],
                 logit:bool = True,
                 anova:bool = False,
                 allow_warnings:bool = False,
                 sequence_context:int = 0,
                 sequence_context_weights:str = "uniform",
                 min_coverage:int = 30,
                 min_ref_length:int = 100,
                 downsample_high_coverage:int = 5000,
                 max_invalid_kmers_freq:float = 0.1,
                 select_ref_id:list = [],
                 exclude_ref_id:list = [],
                 nthreads:int = 3,
                 progress:bool = False):

        """
        Initialise a `SampComp` object and generate a whitelist of references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        Args:
        * input_db_path
            Path to the SQLite database file with event-aligned read/kmer data
        * output_db_path
            Path to the SQLite database file for storing results
        * sample_dict
            Dictionary containing lists of (unique) sample names, grouped by condition
            Example: d = {"control": ["C1", "C2"], "treatment": ["T1", "T2"]}
        * fasta_fn
            Path to a fasta file corresponding to the reference used for read alignment.
            Not needed if 'whitelist' argument is provided.
        * overwrite
            If the output database already exists, overwrite it with a new database?
            By default, new data will be added to previous data.
        * whitelist
            Whitelist object previously generated with nanocompore Whitelist.
            If not given, will be automatically generated.
        * comparison_methods
            Statistical method to compare the two samples (mann_whitney or MW, kolmogorov_smirnov or KS, student_t or ST, gaussian_mixture_model or GMM).
            This can be a list or a comma separated string. {MW,KS,ST,GMM}
        * logit
            Force logistic regression even if we have less than 2 replicates in any condition.
        * allow_warnings
            If True runtime warnings during the ANOVA tests don't raise an error.
        * sequence_context
            Extend statistical analysis to contigous adjacent base if available.
        * sequence_context_weights
            type of weights to used for combining p-values. {uniform,harmonic}
        * min_coverage
            minimal read coverage required in all sample.
        * min_ref_length
            minimal length of a reference transcript to be considered in the analysis
        * downsample_high_coverage
            For reference with higher coverage, downsample by randomly selecting reads.
        * max_invalid_kmers_freq
            maximum frequency of NNNNN, mismatching and missing kmers in reads.
        * select_ref_id
            if given, only reference ids in the list will be selected for the analysis.
        * exclude_ref_id
            if given, refid in the list will be excluded from the analysis.
        * nthreads
            Number of threads (two are used for reading and writing, all the others for parallel processing).
        * progress
            Display a progress bar during execution
        """
        logger.info("Checking and initialising SampComp")

        # Save init options in dict for later
        log_init_state(loc=locals())

        # Check eventalign_dict file paths and labels
        check_sample_dict(sample_dict)
        logger.debug(sample_dict)

        # Check threads number
        if nthreads < 3:
            raise NanocomporeError("The minimum number of threads is 3")

        # Parse comparison methods
        if comparison_methods:
            if type(comparison_methods) == str:
                comparison_methods = comparison_methods.split(",")
            for i, method in enumerate(comparison_methods):
                method = method.upper()
                if method in ["MANN_WHITNEY", "MW"]:
                    comparison_methods[i] = "MW"
                elif method in ["KOLMOGOROV_SMIRNOV", "KS"]:
                    comparison_methods[i] = "KS"
                elif method in ["STUDENT_T", "ST"]:
                    comparison_methods[i] = "ST"
                elif method in ["GAUSSIAN_MIXTURE_MODEL", "GMM"]:
                    comparison_methods[i] = "GMM"
                else:
                    raise NanocomporeError("Invalid comparison method {}".format(method))

        if not whitelist:
            whitelist = Whitelist(input_db_path,
                                  sample_dict,
                                  fasta_fn,
                                  min_coverage = min_coverage,
                                  min_ref_length = min_ref_length,
                                  downsample_high_coverage = downsample_high_coverage,
                                  max_invalid_kmers_freq = max_invalid_kmers_freq,
                                  select_ref_id = select_ref_id,
                                  exclude_ref_id = exclude_ref_id)
        elif not isinstance(whitelist, Whitelist):
            raise NanocomporeError("Whitelist is not valid")

        self.__output_db_path = output_db_path
        db_create_mode = DBCreateMode.OVERWRITE if overwrite else DBCreateMode.CREATE_MAYBE
        with DataStore_SampComp(self.__output_db_path, db_create_mode) as db:
            db.store_whitelist(whitelist)
        # TODO: move this to '__call__'?

        # Set private args from whitelist args
        self.__min_coverage = whitelist._Whitelist__min_coverage
        self.__downsample_high_coverage = whitelist._Whitelist__downsample_high_coverage
        self.__max_invalid_kmers_freq = whitelist._Whitelist__max_invalid_kmers_freq

        # Save private args
        self.__input_db_path = input_db_path
        self.__sample_dict = sample_dict
        self.__fasta_fn = fasta_fn
        self.__whitelist = whitelist
        self.__comparison_methods = comparison_methods
        self.__logit = logit
        self.__anova = anova
        self.__allow_warnings = allow_warnings
        self.__sequence_context = sequence_context
        self.__sequence_context_weights = sequence_context_weights
        self.__nthreads = nthreads - 2
        self.__progress = progress

        # Get number of samples
        self.__n_samples = 0
        for samples in sample_dict.values():
            self.__n_samples += len(samples)


    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Starting data processing")
        # Init Multiprocessing variables
        in_q = mp.Queue(maxsize = 100)
        out_q = mp.Queue(maxsize = 100)
        error_q = mp.Queue()

        # Define processes
        ps_list = []
        ps_list.append(mp.Process(target=self.__list_refid, args=(in_q, error_q)))
        for i in range(self.__nthreads):
            ps_list.append(mp.Process(target=self.__process_references, args=(in_q, out_q, error_q)))
        ps_list.append(mp.Process(target=self.__write_output_to_db, args=(out_q, error_q)))

        # Start processes and monitor error queue
        try:
            # Start all processes
            for ps in ps_list:
                ps.start()

            # Monitor error queue
            for tb in iter(error_q.get, None):
                logger.trace("Error caught from error_q")
                raise NanocomporeError(tb)

            # Soft processes and queues stopping
            for ps in ps_list:
                ps.join()
            for q in (in_q, out_q, error_q):
                q.close()

        # Catch error, kill all processed and reraise error
        except Exception as E:
            logger.error("An error occured. Killing all processes and closing queues\n")
            try:
                for ps in ps_list:
                    ps.terminate()
                for q in (in_q, out_q, error_q):
                    q.close()
            except:
                logger.error("An error occured while trying to kill processes\n")
            raise E


    def process_transcript(self, tx_id, whitelist_reads):
        """Process a transcript given filtered reads from Whitelist"""
        logger.debug(f"Processing transcript: {tx_id}")

        # Kmer data from whitelisted reads from all samples for this transcript
        # Structure: kmer position -> condition -> sample -> data
        kmer_data = defaultdict(lambda: {condition:
                                    defaultdict(lambda: {"intensity": [],
                                                    "dwell": [],
                                                    "coverage": 0,
                                                    "kmers_stats": {"valid": 0,
                                                                    # "missing": 0, # TODO: needed?
                                                                    "NNNNN": 0,
                                                                    "mismatching": 0}})
                                    for condition in self.__sample_dict})
        n_reads = n_kmers = 0

        # Read kmer data from database
        with DataStore_EventAlign(self.__input_db_path) as db:
            for cond_lab, sample_dict in whitelist_reads.items():
                for sample_id, read_ids in sample_dict.items():
                    if not read_ids: continue # TODO: error?
                    n_reads += len(read_ids)
                    values = ", ".join([str(read_id) for read_id in read_ids])
                    query = f"SELECT * FROM kmers WHERE readid IN ({values})"
                    for row in db.cursor.execute(query):
                        n_kmers += 1
                        pos = row["position"]
                        # TODO: check that kmer seq. agrees with FASTA?
                        data = kmer_data[pos][cond_lab][sample_id]
                        data["intensity"].append(row["median"])
                        data["dwell"].append(row["dwell_time"])
                        data["coverage"] += 1
                        status = row["status"]
                        data["kmers_stats"][status] += 1

        logger.debug(f"Data loaded for transcript: {tx_id}")
        test_results = {}
        if self.__comparison_methods:
            random_state = np.random.RandomState(seed=42)
            test_results = txCompare(tx_id,
                                     kmer_data,
                                     random_state=random_state,
                                     methods=self.__comparison_methods,
                                     sequence_context=self.__sequence_context,
                                     sequence_context_weights=self.__sequence_context_weights,
                                     min_coverage= self.__min_coverage,
                                     allow_warnings=self.__allow_warnings,
                                     logit=self.__logit,
                                     anova=self.__anova)

        # Remove 'default_factory' functions from 'kmer_data' to enable pickle/multiprocessing
        kmer_data.default_factory = None
        for pos_dict in kmer_data.values():
            for cond_dict in pos_dict.values():
                cond_dict.default_factory = None

        return {"kmer_data": kmer_data, "test_results": test_results,
                "n_reads": n_reads, "n_kmers": n_kmers}


    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#
    def __list_refid(self, in_q, error_q):
        """Add valid refid from whitelist to input queue to dispatch the data among the workers"""
        n_tx = 0
        try:
            for ref_id, ref_dict in self.__whitelist:
                logger.debug("Adding {} to in_q".format(ref_id))
                in_q.put((ref_id, ref_dict))
                n_tx += 1

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in Reader")
            error_q.put(traceback.format_exc())

        # Deal poison pills
        finally:
            for i in range (self.__nthreads):
                in_q.put(None)
            logger.debug("Parsed transcripts:{}".format(n_tx))


    def __process_references(self, in_q, out_q, error_q):
        """
        Consume ref_id, agregate intensity and dwell time at position level and
        perform statistical analyses to find significantly different regions
        """
        n_tx = n_reads = n_kmers = 0
        try:
            logger.debug("Worker thread started")
            # Process references in input queue
            for ref_id, ref_dict in iter(in_q.get, None):
                logger.debug(f"Worker thread processing new item from in_q: {ref_id}")
                results = self.process_transcript(ref_id, ref_dict)
                n_tx += 1
                n_reads += results["n_reads"]
                n_kmers += results["n_kmers"]

                # Add the current read details to queue
                logger.debug(f"Adding '{ref_id}' to out_q")
                out_q.put((ref_id, results["kmer_data"], results["test_results"]))

        # Manage exceptions and add error trackback to error queue
        except Exception as e:
            logger.error("Error in Worker")
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pill and close file pointer
        finally:
            logger.debug(f"Processed {n_tx} transcripts, {n_reads} reads, {n_kmers} kmers")
            logger.debug("Adding poison pill to out_q")
            out_q.put(None)


    def __write_output_to_db(self, out_q, error_q):
        n_tx = 0
        try:
            # Database was already created earlier to store the whitelist!
            with DataStore_SampComp(self.__output_db_path, DBCreateMode.MUST_EXIST) as db:
                # Iterate over the counter queue and process items until all poison pills are found
                for _ in range(self.__nthreads):
                    for ref_id, kmer_data, test_results in iter(out_q.get, None):
                        logger.debug("Writer thread storing transcript %s" % ref_id)
                        db.store_test_results(ref_id, test_results)
                        n_tx += 1
        except Exception:
            logger.error("Error in writer thread")
            error_q.put(traceback.format_exc())
        finally:
            logger.info(f"All done. Transcripts processed: {n_tx}")
            # Kill error queue with poison pill
            error_q.put(None)


    def __write_output(self, out_q, error_q):
        # Get results out of the out queue and write in shelve
        pvalue_tests = set()
        ref_id_list = []
        n_tx = n_pos = 0
        try:
            with shelve.open(self.__db_fn, flag='n') as db, \
                 tqdm(total=len(self.__whitelist), unit=" Processed References",
                      disable=not self.__progress) as pbar:
                # Iterate over the counter queue and process items until all poison pills are found
                for _ in range(self.__nthreads):
                    for ref_id, kmer_data, test_results in iter(out_q.get, None):
                        ref_id_list.append(ref_id)
                        logger.debug("Writer thread writing %s"%ref_id)
                        # Get pvalue fields available in analysed data
                        for res_dict in test_results.values():
                            for res in res_dict.keys():
                                if "pvalue" in res:
                                    n_pos += 1
                                    pvalue_tests.add(res)
                        # Write results in a shelve db
                        db[ref_id] = (kmer_data, test_results)
                        pbar.update(1)
                        n_tx += 1

                # Write list of refid
                db["__ref_id_list"] = ref_id_list

                # Write metadata
                db["__metadata"] = {
                    "package_name": pkg.__version__,
                    "package_version": pkg.__name__,
                    "timestamp": str(datetime.datetime.now()),
                    "comparison_methods": self.__comparison_methods,
                    "pvalue_tests": sorted(list(pvalue_tests)),
                    "sequence_context": self.__sequence_context,
                    "min_coverage": self.__min_coverage,
                    "n_samples": self.__n_samples}

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.error("Error in Writer")
            error_q.put(traceback.format_exc())

        finally:
            logger.debug(f"Wrote {n_tx} transcripts, {n_pos} valid positions")
            logger.info(f"All done. Transcripts processed: {n_tx}")
            # Kill error queue with poison pill
            error_q.put(None)
