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
import numpy as np
from pyfaidx import Fasta
from statsmodels.stats.multitest import multipletests

# Local package
from nanocompore.common import *
from nanocompore.DataStore import *
from nanocompore.Whitelist import Whitelist
from nanocompore.TxComp import TxComp
from nanocompore.SampCompDB import SampCompDB

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
                 fasta_fn:str,
                 master_db:str = "nanocompore.db",
                 univariate_test:str = "KS", # or: "MW", "ST"
                 fit_gmm:bool = True,
                 gmm_test:str = "logit", # or: "anova"
                 store_gmm_components:str = "none",
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
                 nthreads:int = 3):
        """
        Initialise a `SampComp` object and generate a whitelist of references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        Args:
        * input_dir
            Path to the directory containing input data ("eventalign_collapse" output)
        * fasta_fn
            Path to a fasta file corresponding to the reference used for read alignment.
        * master_db
            Filename of the master database
        * univariate_test
            Statistical test to compare the two conditions ('MW' for Mann-Whitney, 'KS' for Kolmogorov-Smirnov or 'ST' for Student's t), or empty for no test.
        * fit_gmm
            Fit a Gaussian mixture model (GMM) to the intensity/dwell-time distribution?
        * gmm_test
            Method to compare samples based on the GMM ('logit' or 'anova'), or empty for no comparison.
        * store_gmm_components
            Store parameters of any fitted GMMs in the transcript databases? ('none', 'best' or 'all')
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
            For transcripts with high coverage, downsample by selecting reads with most valid kmers.
        * max_invalid_kmers_freq
            maximum frequency of NNNNN, mismatching and missing kmers in reads.
        * significance_thresholds
            Dictionary of significance thresholds used for filtering results by p-value. Only results that meet at least one of the thresholds will be written to the master database (and from there exported). If empty, no filtering is done.
        * select_ref_id - TODO: implement
            if given, only reference ids in the list will be selected for the analysis.
        * exclude_ref_id - TODO: implement
            if given, refid in the list will be excluded from the analysis.
        * nthreads
            Number of threads (two are used for reading and writing, all the others for parallel processing).
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
        if store_gmm_components not in ["none", "best", "all"]:
            raise NanocomporeError(f"Invalid choice for 'store_gmm_components': {store_gmm_components}")

        # Test if Fasta can be opened
        try:
            with Fasta(fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        self._min_transcript_length = min_transcript_length
        self._min_coverage = min_coverage
        self._max_invalid_kmers_freq = max_invalid_kmers_freq
        self._downsample_high_coverage = downsample_high_coverage

        # Prepare database query once
        # TODO: move this to 'DataStore_transcript'?
        subquery = "SELECT id AS reads_id, sampleid, dwelltime_log10_mean, dwelltime_log10_sd FROM reads WHERE pass_filter = 1"
        columns = "sampleid, position, intensity, dwelltime_log10, dwelltime_log10_mean, dwelltime_log10_sd"
        # select only valid kmers (status 0):
        self._kmer_query = f"SELECT {columns} FROM kmers INNER JOIN ({subquery}) ON readid = reads_id WHERE statusid = 0"

        # Save private args
        self._input_dir = input_dir
        self._master_db_path = os.path.join(input_dir, master_db)
        self._nthreads = nthreads - 2 # subtract main thread and writer thread

        # parameters only needed for TxComp:
        self._txcomp_params = {"sequence_context": sequence_context,
                               "sequence_context_weighting": sequence_context_weighting,
                               "allow_anova_warnings": allow_anova_warnings}

        # used to update databases and to adjust p-values:
        self._univariate_test = univariate_test
        self._fit_gmm = fit_gmm
        self._gmm_test = gmm_test if fit_gmm else ""
        self._store_gmm_components = store_gmm_components if fit_gmm else "none"
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
                                    self._max_invalid_kmers_freq,
                                    downsample_high_coverage=self._downsample_high_coverage)

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
        in_q = mp.Queue()
        out_q = mp.Queue(maxsize = 100)
        error_q = mp.Queue()
        # Define processes
        ps_list = [mp.Process(target=self.__process_transcripts, args=(in_q, out_q, error_q))
                   for _ in range(self._nthreads)]
        ps_list.append(mp.Process(target=self.__write_output_to_db, args=(out_q, error_q)))

        logger.debug("Querying transcripts and storing parameters in master DB")
        with DataStore_master(self._master_db_path, DBCreateMode.MUST_EXIST) as db:
            # write parameters to DB (special handling of 'significance_thresholds' dict):
            thresholds = {}
            for k, v in self._significance_thresholds.items():
                thresholds[f"significance_threshold:{k}"] = v
            db.store_parameters("SC", fasta_fn=self._fasta_fn, univariate_test=self._univariate_test,
                                fit_gmm=self._fit_gmm, gmm_test=self._gmm_test,
                                sequence_context=self._txcomp_params["sequence_context"],
                                sequence_context_weighting=self._txcomp_params["sequence_context_weighting"],
                                min_coverage=self._min_coverage, min_transcript_length=self._min_transcript_length,
                                downsample_high_coverage=self._downsample_high_coverage,
                                max_invalid_kmers_freq=self._max_invalid_kmers_freq, **thresholds)
            # enqueue transcripts for processing:
            # DB query needs to finish before processing starts, otherwise master DB will be locked to updates!
            for row in db.cursor.execute("SELECT id, name, subdir FROM transcripts"):
                in_q.put((row["id"], row["name"], row["subdir"]))
        logger.info(f"{in_q.qsize()} transcripts scheduled for processing")
        # Deal poison pills to signal end of input
        for i in range(self._nthreads):
            in_q.put(None)

        # Start processes and monitor error queue
        try:
            # Start all processes
            for ps in ps_list:
                ps.start()
            # Monitor error queue
            for trace in iter(error_q.get, None):
                logger.error("Error caught from error_q")
                raise NanocomporeError(trace)
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
        if self._univariate_test or self._gmm_test:
            logger.info("Running multiple testing correction")
            self.__adjust_pvalues()
            # context-based p-values are not independent tests, so adjust them separately:
            if self._sequence_context:
                self.__adjust_pvalues(sequence_context=True)


    def process_transcript(self, tx_name, subdir):
        """Process data from one transcript"""
        logger.debug(f"Processing transcript: {tx_name}")

        db_path = os.path.join(self._input_dir, subdir, tx_name + ".db")
        if not os.path.exists(db_path):
            logger.error(f"Transcript database not found: {db_path}")
            # TODO: exception?
            return (3, None)
        tx_ok = self._whitelist(db_path)
        if not tx_ok:
            logger.debug("Coverage too low - skipping this transcript")
            return (2, None)

        # Kmer data from whitelisted reads from all samples for this transcript
        # Structure: kmer position -> sample -> data
        kmer_data = defaultdict(lambda: {sample_id: {"intensity": [],
                                                "dwelltime": [],
                                                "coverage": 0}
                                    for sample_id in self._condition_lookup})

        n_kmers = 0 # TODO: count number of reads? (adds overhead)
        # Read kmer data from database
        with DataStore_transcript(db_path, tx_name, subdir) as db:
            # Get transcript-wide intensity stats for scaling:
            db.cursor.execute("SELECT avg(intensity_mean), avg(intensity_sd) FROM reads WHERE pass_filter = 1")
            intensity_mean, intensity_sd = db.cursor.fetchone()
            for row in db.cursor.execute(self._kmer_query):
                n_kmers += 1
                sample_id = row["sampleid"]
                pos = row["position"]
                data = kmer_data[pos][sample_id]
                # Standardize intensity and (log10) dwell time to z-scores:
                data["intensity"].append(self.standard_scale(row["intensity"],
                                                             intensity_mean, intensity_sd))
                data["dwelltime"].append(self.standard_scale(row["dwelltime_log10"],
                                                             row["dwelltime_log10_mean"], row["dwelltime_log10_sd"]))
                data["coverage"] += 1

        logger.trace(f"Data loaded for transcript: {tx_name}")
        if self._tx_compare:
            test_results, n_univariate_tests, n_gmm_tests = self._tx_compare(tx_name, kmer_data)
            # Write complete results to transcript database:
            logger.trace("Writing test results to database")
            with DataStore_transcript(db_path, tx_name, subdir) as db:
                db.create_stats_tables(self._fit_gmm, self._store_gmm_components, self._sequence_context)
                db.store_test_results(test_results)
            if self._significance_thresholds: # keep only significant results for "master" DB
                self.__filter_test_results(test_results)
        else:
            test_results = n_univariate_tests = n_gmm_tests = None

        results = {"test_results": test_results, "n_kmers": n_kmers,
                   "n_univariate_tests": n_univariate_tests, "n_gmm_tests": n_gmm_tests}
        return (0, results)


    @staticmethod
    def standard_scale(x, mean, sd):
        return (x - mean) / sd


    def __filter_test_results(self, test_results):
        logger.trace("Filtering test results for significance")
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


    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHODS~~~~~~~~~~~~~~#

    def __process_transcripts(self, in_q, out_q, error_q):
        """
        Consume one transcript, aggregate intensity and dwell time at position level and
        perform statistical analyses to find significantly different regions
        """
        n_tx = n_kmers = 0
        try:
            logger.debug("Worker thread started")
            # Fasta class may not be multiprocessing-safe, so initialise here (in each process):
            with Fasta(self._fasta_fn) as fasta:
                # Process references in input queue
                for tx_id, tx_name, subdir in iter(in_q.get, None):
                    n_tx += 1
                    logger.trace(f"Worker thread processing new item from in_q: {tx_name}")
                    if len(fasta[tx_name]) < self._min_transcript_length:
                        logger.debug(f"Skipping short transcript '{tx_name}'")
                        out_q.put((tx_id, tx_name, 1, None))
                    else:
                        status, results = self.process_transcript(tx_name, subdir)
                        if results:
                            n_kmers += results["n_kmers"]
                        out_q.put((tx_id, tx_name, status, results))
                    logger.trace(f"Added '{tx_name}' to out_q")

        # Manage exceptions and add error traceback to error queue
        except Exception as e:
            logger.error("Error in Worker")
            error_q.put(NanocomporeError(traceback.format_exc()))

        # Deal poison pill
        finally:
            logger.debug(f"Worker thread processed {n_tx} transcripts, {n_kmers} kmers")
            logger.trace("Adding poison pill to out_q")
            out_q.put(None)


    def __write_output_to_db(self, out_q, error_q):
        n_tx = 0
        try:
            with DataStore_master(self._master_db_path) as db:
                db.init_test_results(bool(self._univariate_test), bool(self._gmm_test), self._sequence_context)
                # Iterate over the output queue and process items until all poison pills are found
                for _ in range(self._nthreads):
                    for tx_id, tx_name, status, results in iter(out_q.get, None):
                        logger.trace(f"Writer thread storing transcript: {tx_name}")
                        db.store_test_results(tx_id, status, results)
                        n_tx += 1
        except Exception:
            logger.error("Error in writer thread")
            error_q.put(traceback.format_exc())
        finally:
            logger.info(f"All done. Transcripts processed: {n_tx}")
            # Kill error queue with poison pill
            error_q.put(None)


    # TODO: move this to 'DataStore_SampComp'?
    def __adjust_pvalues(self, sequence_context=False):
        """Perform multiple testing correction of p-values and update database"""
        if not self._univariate_test and not self._gmm_test:
            return
        pval_columns = ["rowid"] # TODO: use explicit "id" column once defined
        ntest_columns = []
        # for "context-averaged" p-values, add a suffix to the column names:
        col_suffix = "_context" if sequence_context else ""
        if self._univariate_test:
            pval_columns += [f"intensity_pvalue{col_suffix}", f"dwelltime_pvalue{col_suffix}"]
            ntest_columns.append("SUM(n_univariate_tests)")
        if self._gmm_test:
            pval_columns.append(f"gmm_pvalue{col_suffix}")
            ntest_columns.append("SUM(n_gmm_tests)")
        with DataStore_master(self._master_db_path) as db:
            # prepare database for updates:
            for col in pval_columns[1:]:
                db.add_or_reset_column("test_results", f"adj_{col}", "REAL")
            # get p-values from database:
            pvalues = []
            index = []
            columns = ", ".join(pval_columns)
            try:
                for row in db.cursor.execute(f"SELECT {columns} FROM test_results"):
                    for pv_col in row.keys()[1:]:
                        pv = row[pv_col]
                        # correction function doesn't handle NaN values, so skip those:
                        if (pv is not None) and not np.isnan(pv):
                            pvalues.append(pv)
                            index.append({"id": row["rowid"], "pv_col": pv_col})
            except:
                logger.error("Error reading p-values from table 'test_results'")
                raise
            logger.debug(f"Number of p-values to correct for multiple testing: {len(pvalues)}")
            if not pvalues:
                return
            # get total number of tests:
            columns = ", ".join(ntest_columns)
            try:
                row = db.cursor.execute(f"SELECT {columns} FROM transcripts").fetchone()
                n_tests = sum(row)
            except:
                logger.error("Error reading number of statistical tests from table 'transcripts'")
            logger.debug(f"Number of statistical tests performed: {n_tests}")
            # perform multiple testing correction:
            adjusted = self.fdr_adjust(pvalues, n_tests)
            assert len(pvalues) == len(adjusted)
            # update database:
            # sqlite module can't handle numpy float64 values, so convert to floats using "tolist":
            for ind, adj_pv in zip(index, adjusted.tolist()):
                query = "UPDATE test_results SET adj_{pv_col} = ? WHERE rowid = {id}".format_map(ind)
                try:
                    db.cursor.execute(query, (adj_pv, ))
                except:
                    logger.error("Error updating adjusted p-value for ID {id}".format_map(ind))
                    raise


    @staticmethod
    def fdr_adjust(pvalues, n_tests=0):
        n_tests = max(n_tests, len(pvalues))
        if (len(pvalues) == 0) or (n_tests == 1):
            return pvalues
        # false discovery rate p-value adjustment (Benjamini-Hochberg method):
        # see https://stackoverflow.com/a/33532498 - adapted for "n_tests"
        # results have been checked against R function 'p.adjust(method="fdr")'
        p = np.asfarray(pvalues)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(n_tests) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]


    def reset_databases(self):
        with DataStore_master(self._master_db_path) as master:
            master.reset_SampComp()
            for row in master.cursor.execute("SELECT id, name, subdir FROM transcripts"):
                db_path = os.path.join(self._input_dir, row["subdir"], row["name"] + ".db")
                with DataStore_transcript(db_path, row["name"], row["id"]) as db:
                    db.reset_SampComp()
        logger.info("All databases reset to pre-SampComp state")
