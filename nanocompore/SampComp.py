# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Disable multithreading for MKL and openBlas
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# Std lib
import logging
from collections import OrderedDict, Counter
import shelve
import multiprocessing as mp
from warnings import warn

# Third party
from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta

# Local package
from nanocompore.common import counter_to_str, access_file, NanocomporeError, NanocomporeWarning, numeric_cast_list
from nanocompore.Whitelist import Whitelist
from nanocompore.TxComp import txCompare
from nanocompore.SampCompDB import SampCompDB

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampComp (object):
    """ Init analysis and check args"""

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

    def __init__(self,
        eventalign_fn_dict,
        output_db_fn,
        fasta_fn,
        bed_fn=None,
        whitelist=None,
        comparison_method = None,
        sequence_context = 0,
        min_coverage = 10,
        downsample_high_coverage = None,
        max_invalid_kmers_freq = 0.1,
        select_ref_id = [],
        exclude_ref_id = [],
        nthreads = 4,
        logLevel = "info"):

        """
        eventalign_fn_dict: Multilevel dictionnary indicating the condition_label, sample_label and file name of the eventalign_collapse output
            example d = {"S1": {"R1":"path1.tsv", "R2":"path2.tsv"}, "S2": {"R1":"path3.tsv", "R2":"path4.tsv"}}
            2 conditions are expected, and at least 2 sample replicates are highly recomended per condition
        output_db_fn: Path where to write the result database
        fasta_fn: Path to a fasta file corresponding to the reference used for read alignemnt
        bed_fn: Path to a BED file containing the annotation of the transcriptome used as reference when mapping
        whitelist: Whitelist object previously generated with nanocompore Whitelist. If not given, will be automatically generated
        comparison_method: Statistical method to compare the 2 samples (mann_whitney, kolmogorov_smirnov, t_test, gmm).
            This can be a list or a comma separated string
        sequence_context: Extend statistical analysis to contigous adjacent base if available
        min_cov: minimal coverage required in all sample
        downsample_high_coverage: For reference with higher coverage, downsample by randomly selecting reads.
        max_invalid_kmers_freq: maximum frequency of NNNNN, mismatching and missing kmers in reads
        select_ref_id: if given, only reference ids in the list will be selected for the analysis
        exclude_ref_id: if given, refid in the list will be excluded from the analysis
        nthreads: Number of threads (two are used for reading and writing, all the others for processing in parallel).
        logLevel: Set the log level. Valid values: warning, info, debug
        """
        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise SampComp and checks options")

        # Check that the number of condition is 2 and raise a warning if there are less than 2 replicates per conditions
        if len(eventalign_fn_dict) != 2:
            raise NanocomporeError("2 conditions are expected. Found {}".format(len(eventalign_fn_dict)))
        for cond_lab, sample_dict in eventalign_fn_dict.items():
            if len(sample_dict) == 1:
                warn (NanocomporeWarning ("Only 1 replicate found for condition {}. This is not recomended".format(cond_lab)))

        # Check args
        for sample_dict in eventalign_fn_dict.values():
            for fn in sample_dict.values():
                if not access_file (fn):
                    raise NanocomporeError("Cannot access eventalign_collapse file {}".format(idx_fn))

        if nthreads < 3:
            raise NanocomporeError("Number of threads not valid")

        # Parse comparison methods
        if comparison_method:
            if type (comparison_method) == str:
                comparison_method = comparison_method.split(",")
            for i, method in enumerate (comparison_method):
                method = method.upper()
                if method in ["MANN_WHITNEY", "MW"]:
                    comparison_method[i]="MW"
                elif method in ["KOLMOGOROV_SMIRNOV", "KS"]:
                    comparison_method[i]="KS"
                elif method in ["T_TEST", "TT"]:
                    comparison_method[i]="TT"
                elif method in ["GAUSSIAN_MIXTURE_MODEL", "GMM"]:
                    comparison_method[i]="GMM"
                else:
                    raise NanocomporeError("Invalid comparison method {}".format(method))

        if not whitelist:
            whitelist = Whitelist (
                eventalign_fn_dict = eventalign_fn_dict,
                fasta_fn = fasta_fn,
                min_coverage = min_coverage,
                downsample_high_coverage = downsample_high_coverage,
                max_invalid_kmers_freq = max_invalid_kmers_freq,
                select_ref_id = select_ref_id,
                exclude_ref_id = exclude_ref_id,
                logLevel = logLevel)
        elif not isinstance (whitelist, Whitelist):
            raise NanocomporeError("Whitelist is not valid")

        # Set private args from whitelist args
        self.__min_coverage = whitelist._Whitelist__min_coverage
        self.__downsample_high_coverage = whitelist._Whitelist__downsample_high_coverage
        self.__max_invalid_kmers_freq = whitelist._Whitelist__max_invalid_kmers_freq

        # Save private args
        self.__eventalign_fn_dict = eventalign_fn_dict
        self.__output_db_fn = output_db_fn
        self.__fasta_fn = fasta_fn
        self.__bed_fn = bed_fn
        self.__whitelist = whitelist
        self.__comparison_methods = comparison_method
        self.__sequence_context = sequence_context
        self.__nthreads = nthreads - 2
        self.__logLevel = logLevel

        # Get number of samples
        n = 0
        for sample_dict in self.__eventalign_fn_dict.values():
            for sample_lab in sample_dict.keys():
                n+=1
        self.__n_samples = n

    def __call__ (self):
        """Run analysis"""

        logger.info("Start data processing")
        # Init Multiprocessing variables
        in_q = mp.Queue(maxsize = 100)
        out_q = mp.Queue(maxsize = 100)

        # Define processes
        ps_list = []
        ps_list.append(mp.Process(target=self.__list_refid, args=(in_q,)))
        for i in range(self.__nthreads):
            ps_list.append(mp.Process(target=self.__process_references, args=(in_q, out_q)))
        ps_list.append(mp.Process(target=self.__write_output, args=(out_q,)))

        # Start processes and block until done
        try:
            for ps in ps_list:
                ps.start()
            for ps in ps_list:
                ps.join()

        # Kill processes if early stop
        except (BrokenPipeError, KeyboardInterrupt) as E:
            if self.verbose: stderr_print("Early stop. Kill processes\n")
            for ps in ps_list:
                ps.terminate()

        # Return database wrapper object
        return SampCompDB(db_fn=self.__output_db_fn, fasta_fn=self.__fasta_fn, bed_fn=self.__bed_fn)

    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#
    def __list_refid (self, in_q):
        # Add refid to inqueue to dispatch the data among the workers
        for ref_id, ref_dict in self.__whitelist:
            logger.debug("Adding {} to in_q".format(ref_id))
            in_q.put((ref_id, ref_dict))

        # Add 1 poison pill for each worker thread
        for i in range(self.__nthreads):
            logger.debug("Adding poison pill to in_q")
            in_q.put(None)

    def __process_references (self, in_q, out_q):
        """
        Consume ref_id, agregate intensity and dwell time at position level and
        perform statistical analyses to find significantly different regions
        """
        logger.debug("Worker thread started")

        try:
            # Open all files for reading. File pointer are stored in a dict matching the ref_dict entries
            fp_dict = self.__eventalign_fn_open()

            # Process refid in input queue
            for ref_id, ref_dict in iter(in_q.get, None):
                logger.debug("Worker thread processing new item from in_q: {}".format(ref_id))

                # Create an empty dict for all positions first
                ref_pos_dict = self.__make_ref_pos_dict (ref_id)

                for cond_lab, sample_dict in ref_dict["data"].items():
                    for sample_lab, read_list in sample_dict.items():
                        fp = fp_dict[cond_lab][sample_lab]

                        for read in read_list:

                            # Move to read, save read data chunk and reset file pointer
                            fp.seek (read["byte_offset"])
                            line_list = fp.read (read["byte_len"]).split("\n")
                            fp.seek (0)

                            # Check read_id ref_id concordance between index and data file
                            header = numeric_cast_list(line_list[0][1:].split("\t"))
                            if not header[0] == read["read_id"] or not header[1] == read["ref_id"]:
                                raise NanocomporeError ("Index and data files are not matching")

                            # Extract col names from second line
                            col_names = line_list[1].split("\t")

                            # Parse data files kmers per kmers
                            prev_pos = None
                            for line in line_list[2:]:
                                # Transform line to dict and cast str numbers to actual numbers
                                kmer = dict (zip (col_names, numeric_cast_list (line.split("\t"))))

                                # Check if ref position is in the whitelist valid intervals and add kmer data if so
                                pos = kmer["ref_pos"]

                                # Fill in the missing positions
                                if prev_pos and pos-prev_pos > 1:
                                    for missing_pos in range(prev_pos+1, pos):
                                        ref_pos_dict[missing_pos]["data"][cond_lab][sample_lab]["events_stats"]["missing"] += 1

                                # And then fill dict with the current pos values
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["intensity"].append(kmer["median"])
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["dwell"].append(kmer["n_signals"])
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["coverage"] += 1

                                # Also fill in with normalised position event stats
                                n_valid = (kmer["n_events"]-(kmer["NNNNN_events"]+kmer["mismatching_events"])) / kmer["n_events"]
                                n_NNNNN = kmer["NNNNN_events"] / kmer["n_events"]
                                n_mismatching = kmer["mismatching_events"] / kmer["n_events"]
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["events_stats"]["valid"] += n_valid
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["events_stats"]["NNNNN"] += n_NNNNN
                                ref_pos_dict[pos]["data"][cond_lab][sample_lab]["events_stats"]["mismatching"] += n_mismatching

                                prev_pos = pos

                if self.__comparison_methods:
                    ref_pos_dict=txCompare(
                        data=ref_pos_dict,
                        methods=self.__comparison_methods,
                        sequence_context=self.__sequence_context,
                        min_coverage=self.__min_coverage,
                        logger=logger,
                        ref=ref_id)

                # Add the current read details to queue
                logger.debug("Adding %s to out_q"%(ref_id))
                out_q.put ((ref_id, ref_pos_dict))

        finally:
            # Add a poison pill in queue
            logger.debug("Adding poison pill to out_q")
            out_q.put (None)
            # close all files
            self.__eventalign_fn_close (fp_dict)

    def __write_output (self, out_q):
        # Get results out of the out queue and write in shelve
        try:
            with shelve.open (self.__output_db_fn, flag='n') as db:
                # Iterate over the counter queue and process items until all poison pills are found
                pbar = tqdm (total = len(self.__whitelist), unit=" Processed References", disable=self.__logLevel in ("warning", "debug"))
                for _ in range (self.__nthreads):
                    for ref_id, ref_pos_dict in iter (out_q.get, None):
                        logger.debug("Writer thread writing %s"%ref_id)
                        # Write results in a shelve db
                        db [ref_id] = ref_pos_dict
                        pbar.update ()
                pbar.close()

                db["__metadata"] = {
                    "comparison_method": self.__comparison_methods,
                    "sequence_context": self.__sequence_context,
                    "min_coverage": self.__min_coverage,
                    "n_samples": self.__n_samples}
        except IOError:
            raise NanocomporeError("Error writing to output db")

    #~~~~~~~~~~~~~~PRIVATE HELPER METHODS~~~~~~~~~~~~~~#
    def __eventalign_fn_open (self):
        fp_dict = OrderedDict()
        for cond_lab, sample_dict in self.__eventalign_fn_dict.items():
            fp_dict[cond_lab] = OrderedDict()
            for sample_lab, fn in sample_dict.items():
                fp_dict[cond_lab][sample_lab] = open(fn, "r")
        return fp_dict

    def __eventalign_fn_close (self, fp_dict):
        for sample_dict in fp_dict.values():
            for fp in sample_dict.values():
                fp.close()

    def __make_ref_pos_dict (self, ref_id):

        ref_pos_dict = OrderedDict()
        with Fasta (self.__fasta_fn) as fasta:
            ref_fasta = fasta [ref_id]
            ref_len = len(ref_fasta)
            ref_seq = str(ref_fasta)+"#####" # Add padding for last kmers

            for pos in range(ref_len):
                pos_dict = OrderedDict()
                pos_dict["ref_kmer"] = ref_seq[pos:pos+5]
                pos_dict["data"] = OrderedDict()
                for cond_lab, s_dict in self.__eventalign_fn_dict.items():
                    pos_dict["data"][cond_lab] = OrderedDict()
                    for sample_lab in s_dict.keys():
                        pos_dict["data"][cond_lab][sample_lab] = {
                            "intensity":[],
                            "dwell":[],
                            "coverage":0,
                            "events_stats":{"missing":0,"valid":0,"NNNNN":0,"mismatching":0}}
                ref_pos_dict[pos] = pos_dict
        return ref_pos_dict
