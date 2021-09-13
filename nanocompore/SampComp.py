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
from statsmodels.stats.multitest import multipletests

# Local package
from nanocompore.common import *
from nanocompore.DataStore import *
from nanocompore.Whitelist import Whitelist
from nanocompore.TxComp import TxComp
from nanocompore.SampCompDB import SampCompDB
import nanocompore as pkg

# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampComp(object):
    """Init analysis and check args"""

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

    # TODO: use enums for univariate and gmm test parameters?
    def __init__(self,
                 input_dir:str,
                 master_db:str = "eventalign_collapse.db",
                 fasta_fn:str = "",
                 univariate_test:str = "KS", # or: "MW", "ST"
                 fit_gmm:bool = True,
                 gmm_test:str = "logit", # or: "anova"
                 allow_anova_warnings:bool = False,
                 sequence_context:int = 0,
                 sequence_context_weighting:str = "uniform",
                 min_coverage:int = 30,
                 min_transcript_length:int = 100,
                 downsample_high_coverage:int = 5000,
                 max_invalid_kmers_freq:float = 0.1,
                 significance_thresholds = {"gmm_pvalue": 0.01},
                 # select_ref_id:list = [],
                 # exclude_ref_id:list = [],
                 nthreads:int = 3,
                 progress:bool = False):

        """
        Initialise a `SampComp` object and generate a whitelist of references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        Args:
        * input_dir
            Path to the directory containing input data ("eventalign_collapse" output)
        * master_db
            Filename of the master database
        * fasta_fn
            Path to a fasta file corresponding to the reference used for read alignment.
            Not needed if 'whitelist' argument is provided.
        * univariate_test
            Statistical test to compare the two conditions ('MW' for Mann-Whitney, 'KS' for Kolmogorov-Smirnov or 'ST' for Student's t), or empty for no test.
        * fit_gmm
            Fit a Gaussian mixture model (GMM) to the intensity/dwell-time distribution?
        * gmm_test
            Method to compare samples based on the GMM ('logit' or 'anova'), or empty for no comparison.
        * allow_anova_warnings
            If True runtime warnings during the ANOVA tests don't raise an error.
        * sequence_context
            Extend statistical analysis to contiguous adjacent bases if available.
        * sequence_context_weighting
            type of weighting to used for combining p-values. {uniform,harmonic}
        * min_coverage
            minimal read coverage required in all samples.
        * min_transcript_length
            minimal length of a reference transcript to be considered in the analysis
        * downsample_high_coverage
            For reference with higher coverage, downsample by randomly selecting reads.
        * max_invalid_kmers_freq
            maximum frequency of NNNNN, mismatching and missing kmers in reads.
        * select_ref_id - TODO: implement
            if given, only reference ids in the list will be selected for the analysis.
        * exclude_ref_id - TODO: implement
            if given, refid in the list will be excluded from the analysis.
        * nthreads
            Number of threads (two are used for reading and writing, all the others for parallel processing).
        * progress
            Display a progress bar during execution
        """
        logger.info("Checking and initialising SampComp")

        # Save init options in dict for later
        log_init_state(loc=locals())

        # Check threads number
        if nthreads < 3:
            raise NanocomporeError("The minimum number of threads is 3")

        # Parse comparison methods
        if univariate_test and (univariate_test not in ["MW", "KS", "ST"]):
            raise NanocomporeError(f"Invalid univariate test {univariate_test}")
        if fit_gmm and gmm_test and (gmm_test not in ["logit", "anova"]):
            raise NanocomporeError(f"Invalid GMM-based test {gmm_test}")

        # Test if Fasta can be opened
        try:
            with Fasta(fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        self._min_transcript_length = min_transcript_length
        self._min_coverage = min_coverage
        self._max_invalid_kmers_freq = max_invalid_kmers_freq

        # Prepare database query once
        # TODO: move this to 'DataStore_transcript'?
        subquery = "SELECT id AS reads_id, sampleid FROM reads WHERE pass_filter = 1"
        if downsample_high_coverage: # choose reads with most valid kmers
            subquery += f" ORDER BY valid_kmers DESC LIMIT {downsample_high_coverage}"
        columns = "sampleid, position, sequenceid, statusid, dwell_time, median"
        # select only valid kmers (status 0):
        self._kmer_query = f"SELECT {columns} FROM kmers INNER JOIN ({subquery}) ON readid = reads_id WHERE statusid = 0"

        # Save private args
        self._input_dir = input_dir
        self._master_db_path = os.path.join(input_dir, master_db)
        self._nthreads = nthreads - 2
        self._progress = progress

        # parameters only needed for TxComp:
        self._txcomp_params = {"sequence_context": sequence_context,
                               "sequence_context_weighting": sequence_context_weighting,
                               "allow_anova_warnings": allow_anova_warnings}

        # used to update databases and to adjust p-values:
        self._univariate_test = univariate_test
        self._fit_gmm = fit_gmm
        self._gmm_test = gmm_test if fit_gmm else ""
        self._sequence_context = sequence_context > 0

        # Cut-offs for filtering final results:
        self._significance_thresholds = significance_thresholds

        # Get sample IDs from database
        with DataStore_master(self._master_db_path, DBCreateMode.MUST_EXIST) as db:
            self._db_samples = db.get_sample_info()
        if len(self._db_samples) != 2:
            raise NanocomporeError(f"Expected two experimental conditions, found {len(self._db_samples)}: {', '.join(self._db_samples.keys())}")
        # Generate lookup dict. for sample IDs -> conditions
        self._condition_lookup = {}
        for cond, samples in self._db_samples.items():
            for sid, _ in samples:
                self._condition_lookup[sid] = cond

        # Initialise the "Whitelist" (filtering) object:
        self._whitelist = Whitelist(self._db_samples,
                                    self._min_coverage,
                                    self._max_invalid_kmers_freq)

        # If statistical tests are requested, initialise the "TxComp" object:
        if self._univariate_test or self._fit_gmm:
            # Need at least two samples per condition for ANOVA on GMM fit:
            if (self._gmm_test == "anova") and not all([len(samples) > 1 for samples in self._db_samples.values()]):
                logger.warning("Not enough replicates for 'anova' GMM test. Switching to 'logit' test.")
                self._gmm_test = "logit"
            random_state = np.random.RandomState(seed=42)
            self._tx_compare = TxComp(random_state,
                                      self._db_samples,
                                      univariate_test=self._univariate_test,
                                      fit_gmm=self._fit_gmm,
                                      gmm_test=self._gmm_test,
                                      min_coverage=self._min_coverage,
                                      **self._txcomp_params)
        else:
            self._tx_compare = None


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
        ps_list.append(mp.Process(target=self.__list_transcripts_check_length, args=(in_q, error_q)))
        for i in range(self._nthreads):
            ps_list.append(mp.Process(target=self.__process_transcripts, args=(in_q, out_q, error_q)))
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

        # Adjust p-values for multiple testing:
        # if self._univariate_test or self._gmm_test:
        #     logger.info("Running multiple testing correction")
        #     self.__adjust_pvalues()
        #     # context-based p-values are not independent tests, so adjust them separately:
        #     if self._sequence_context:
        #         self.__adjust_pvalues(sequence_context=True)


    def process_transcript(self, tx_name, subdir):
        """Process data from one transcript"""
        logger.debug(f"Processing transcript: {tx_name}")

        db_path = os.path.join(self._input_dir, subdir, tx_name + ".db")
        if not os.path.exists(db_path):
            logger.error(f"Transcript database not found: {db_path}")
            # TODO: exception?
            return None
        tx_ok = self._whitelist(db_path)
        if not tx_ok:
            logger.debug("Coverage too low - skipping this transcript")
            return None

        # Kmer data from whitelisted reads from all samples for this transcript
        # Structure: kmer position -> sample -> data
        kmer_data = defaultdict(lambda: {sample_id: {"intensity": [],
                                                "dwell": [],
                                                "coverage": 0}
                                    for sample_id in self._condition_lookup})

        n_kmers = 0 # TODO: count number of reads? (adds overhead)
        # Read kmer data from database
        with DataStore_transcript(db_path, tx_name, subdir) as db:
            for row in db.cursor.execute(self._kmer_query):
                n_kmers += 1
                sample_id = row["sampleid"]
                pos = row["position"]
                data = kmer_data[pos][sample_id]
                data["intensity"].append(row["median"])
                data["dwell"].append(row["dwell_time"])
                data["coverage"] += 1

        logger.debug(f"Data loaded for transcript: {tx_name}")
        if self._tx_compare:
            test_results, n_univariate_tests, n_gmm_tests = self._tx_compare(tx_name, kmer_data)
            # TODO: check "gmm_anova_failed" state of TxComp object
            # Write complete results to transcript database:
            logger.debug("Writing test results to database")
            with DataStore_transcript(db_path, tx_name, subdir) as db:
                db.create_stats_tables(self._fit_gmm, self._sequence_context)
                db.store_test_results(test_results)
            # Keep only significant results for "master" database:
            self.__filter_test_results(test_results)
        else:
            test_results = n_univariate_tests = n_gmm_tests = None

        return {"test_results": test_results, "n_kmers": n_kmers,
                "n_univariate_tests": n_univariate_tests, "n_gmm_tests": n_gmm_tests}


    def __filter_test_results(self, test_results):
        logger.debug("Filtering test results for significance")
        to_remove = []
        for pos, res in test_results.items():
            # Keep a result if it meets any of the p-value thresholds; otherwise, remove
            for pv_key, threshold in self._significance_thresholds.items():
                pvalue = res.get(pv_key)
                if (pvalue is not None) and (pvalue <= threshold):
                    break
            else: # belongs to the "for" loop!
                to_remove.append(pos)
        for pos in to_remove:
            del test_results[pos]


    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#

    def __list_transcripts_check_length(self, in_q, error_q):
        """Add transcripts to input queue to dispatch the data among the workers"""
        n_tx = 0
        short_transcripts = [] # TODO: record short/skipped transcripts in the database?
        pre_queue = []
        try:
            # DB query needs to finish before further processing, otherwise DB will be locked to updates!
            with DataStore_master(self._master_db_path, DBCreateMode.MUST_EXIST) as db, \
                 Fasta(self._fasta_fn) as fasta:
                for row in db.cursor.execute("SELECT * FROM transcripts"):
                    n_tx += 1
                    accession = row["name"]
                    if len(fasta[accession]) < self._min_transcript_length:
                        logger.debug(f"Skipping short transcript '{accession}'")
                        short_transcripts.append((row["id"], accession))
                    else:
                        logger.debug(f"Transcript '{accession}' will be added to processing queue")
                        pre_queue.append((row["id"], accession, row["subdir"]))
            for item in pre_queue:
                in_q.put(item)
        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in Reader")
            error_q.put(traceback.format_exc())
        # Deal poison pills
        finally:
            for i in range(self._nthreads):
                in_q.put(None)
            logger.debug(f"Found {n_tx} transcripts, skipped {len(short_transcripts)}.")


    def __process_transcripts(self, in_q, out_q, error_q):
        """
        Consume ref_id, agregate intensity and dwell time at position level and
        perform statistical analyses to find significantly different regions
        """
        n_tx = n_kmers = 0
        try:
            logger.debug("Worker thread started")
            # Process references in input queue
            for id, tx_name, subdir in iter(in_q.get, None):
                logger.debug(f"Worker thread processing new item from in_q: {tx_name}")
                results = self.process_transcript(tx_name, subdir)
                if not results:
                    continue
                n_tx += 1
                # n_reads += results["n_reads"]
                n_kmers += results["n_kmers"]
                # Add the current read details to queue
                if results["test_results"]:
                    logger.debug(f"Adding '{tx_name}' to out_q")
                    out_q.put((id, tx_name, results))

        # Manage exceptions and add error traceback to error queue
        except Exception as e:
            logger.error("Error in Worker")
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pill and close file pointer
        finally:
            logger.debug(f"Processed {n_tx} transcripts, {n_kmers} kmers")
            logger.debug("Adding poison pill to out_q")
            out_q.put(None)


    def __write_output_to_db(self, out_q, error_q):
        n_tx = 0
        try:
            with DataStore_master(self._master_db_path) as db:
                db.init_test_results(bool(self._univariate_test), bool(self._gmm_test), self._sequence_context)
                # Iterate over the counter queue and process items until all poison pills are found
                for _ in range(self._nthreads):
                    for id, tx_name, results in iter(out_q.get, None):
                        logger.debug("Writer thread storing transcript %s" % tx_name)
                        db.store_test_results(id, results)
                        n_tx += 1
        except Exception:
            logger.error("Error in writer thread")
            error_q.put(traceback.format_exc())
        finally:
            logger.info(f"All done. Transcripts processed: {n_tx}")
            # Kill error queue with poison pill
            error_q.put(None)


    # TODO: move this to 'DataStore_SampComp'?
    def __adjust_pvalues(self, method="fdr_bh", sequence_context=False):
        """Perform multiple testing correction of p-values and update database"""
        db = DataStore_master(self._master_db_path)
        with db:
            pvalues = []
            index = []
            # for "context-averaged" p-values, add a suffix to the column names:
            col_suffix = "_context" if sequence_context else ""
            if self._univariate_test:
                query = f"SELECT id, intensity_pvalue{col_suffix}, dwell_pvalue{col_suffix} FROM kmer_stats"
                try:
                    for row in db.cursor.execute(query):
                        for pv_col in ["intensity_pvalue", "dwell_pvalue"]:
                            pv_col += col_suffix
                            pv = row[pv_col]
                            # "multipletests" doesn't handle NaN values well, so skip those:
                            if (pv is not None) and not np.isnan(pv):
                                pvalues.append(pv)
                                index.append({"table": "kmer_stats", "id_col": "id",
                                              "id": row["id"], "pv_col": pv_col})
                except:
                    logger.error("Error reading p-values from table 'kmer_stats'")
                    raise
            if self._gmm_test:
                pv_col = "test_pvalue" + col_suffix
                query = f"SELECT kmer_statsid, {pv_col} FROM gmm_stats WHERE {pv_col} IS NOT NULL"
                try:
                    for row in db.cursor.execute(query):
                        pv = row[pv_col]
                        # "multipletests" doesn't handle NaN values well, so skip those:
                        if not np.isnan(pv): # 'None' (NULL) values have been excluded in the query
                            pvalues.append(pv)
                            index.append({"table": "gmm_stats", "id_col": "kmer_statsid",
                                          "id": row["kmer_statsid"], "pv_col": pv_col})
                except:
                    logger.error("Error reading p-values from table 'gmm_stats'")
                    raise
            logger.debug(f"Number of p-values for multiple testing correction: {len(pvalues)}")
            if not pvalues:
                return
            adjusted = multipletests(pvalues, method=method)[1]
            assert len(pvalues) == len(adjusted)
            # sqlite module can't handle numpy float64 values, so convert to floats using "tolist":
            for ind, adj_pv in zip(index, adjusted.tolist()):
                query = "UPDATE {table} SET adj_{pv_col} = ? WHERE {id_col} = {id}".format_map(ind)
                try:
                    db.cursor.execute(query, (adj_pv, ))
                except:
                    logger.error("Error updating adjusted p-value for ID {id} in table '{table}'".format_map(ind))
                    raise
