# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
import logging
from collections import OrderedDict, namedtuple
import shelve
import multiprocessing as mp

# Third party
from tqdm import tqdm

# Local package
#from nanocompore.txCompare import txCompare
from nanocompore.common import counter_to_str, access_file, mytqdm, NanocomporeError
from nanocompore.Whitelist import Whitelist
from nanocompore.TxComp import paired_test

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampComp (object):
    """ Produce useful results. => Thanks Tommaso ! That's a very *useful* comment :P """

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#

    def __init__(self,
        s1_fn,
        s2_fn,
        output_db_fn,
        fasta_index_fn=None,
        whitelist=None,
        padj_threshold = 0.1,
        comparison_method = None,
        sequence_context = 0,
        min_coverage = 10,
        downsample_high_coverage = None, ########################################## Not super happy with that
        max_NNNNN_freq = 0.2,
        max_mismatching_freq = 0.2,
        max_missing_freq = 0.2,
        nthreads = 4,
        logLevel = "info"):

        """
        s1_fn: Path to sample 1 eventalign_collapse data file
        s2_fn: Path to sample 2 eventalign_collapse data file
        output_db_fn: Path where to write the result database
        whitelist: Whitelist object previously generated with nanocompore Whitelist. If not given, will be generated
        padj_threshold: Adjusted p-value threshold for reporting sites.
        comparison_method: Statistical method to compare the 2 samples (kmean, mann_whitney, kolmogorov_smirnov, t_test)
        sequence_context: Extend statistical analysis to contigous adjacent base is available
        nthreads: Number of threads (two are used for reading and writing, all the others for processing in parallel).
        logLevel: Set the log level. Valid values: warning, info, debug
        """
        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise and checks options")

        # Check args
        for fn in (s1_fn, s2_fn):
            if not access_file (fn):
                raise NanocomporeError("Cannot access file {}".format(fn))

        if nthreads < 3:
            raise NanocomporeError("Number of threads not valid")

        if not comparison_method in ["kmean", "mann_whitney", "MW", "kolmogorov_smirnov", "KS","t_test", "TT", None]:
            raise NanocomporeError("Invalid comparison method")

        if whitelist:
            if not isinstance (whitelist, Whitelist):
                raise NanocomporeError("Whitelist is not valid")
        else:
            if not fasta_index_fn:
                raise NanocomporeError("fasta_index_fn is required if no whitelist is provided")

            whitelist = Whitelist (
                s1_index_fn = s1_fn+".idx",
                s2_index_fn = s2_fn+".idx",
                fasta_index_fn = fasta_index_fn,
                min_coverage = min_coverage,
                downsample_high_coverage = downsample_high_coverage,
                max_NNNNN_freq = max_NNNNN_freq,
                max_mismatching_freq = max_mismatching_freq,
                max_missing_freq = max_missing_freq,
                logLevel = logLevel)

        # Save private args
        self.__s1_fn = s1_fn
        self.__s2_fn = s2_fn
        self.__output_db_fn = output_db_fn
        self.__whitelist = whitelist
        self.__padj_threshold = padj_threshold
        self.__comparison_method = comparison_method
        self.__sequence_context = sequence_context
        self.__min_coverage = min_coverage
        self.__downsample_high_coverage = downsample_high_coverage
        self.__max_NNNNN_freq = max_NNNNN_freq
        self.__max_mismatching_freq = max_mismatching_freq
        self.__max_missing_freq = max_missing_freq
        self.__nthreads = nthreads - 2
        self.__logLevel = logLevel

        logger.info ("Start data processing")
        # Init Multiprocessing variables
        in_q = mp.Queue (maxsize = 1000)
        out_q = mp.Queue (maxsize = 1000)

        # Define processes
        ps_list = []
        ps_list.append (mp.Process (target=self.__list_refid, args=(in_q,)))
        for i in range (self.__nthreads):
            ps_list.append (mp.Process (target=self.__process_references, args=(in_q, out_q)))
        ps_list.append (mp.Process (target=self.__write_output, args=(out_q,)))

        # Start processes and block until done
        try:
            for ps in ps_list:
                ps.start ()
            for ps in ps_list:
                ps.join ()

        # Kill processes if early stop
        except (BrokenPipeError, KeyboardInterrupt) as E:
            if self.verbose: stderr_print ("Early stop. Kill processes\n")
            for ps in ps_list:
                ps.terminate ()

    #~~~~~~~~~~~~~~PRIVATE MULTIPROCESSING METHOD~~~~~~~~~~~~~~#
    def __list_refid (self, in_q):
        # Add refid to inqueue to dispatch the data among the workers
        for ref_id in self.__whitelist.ref_id_list:
            in_q.put (ref_id)

        # Add 1 poison pill for each worker thread
        for i in range (self.__nthreads):
            in_q.put (None)

    def __process_references (self, in_q, out_q):
        # Consumme ref_id until empty and perform statiscical analysis

        with open (self.__s1_fn) as s1_fp, open (self.__s2_fn) as s2_fp: # More efficient to open only once the files
            for ref_id in iter (in_q.get, None):

                ref_dict = self.__whitelist[ref_id]
                ref_pos_dict = OrderedDict ()

                for interval_start, interval_end in ref_dict["interval_list"]:
                    for i in range (interval_start, interval_end+1):
                        ref_pos_dict[i] = {
                            "S1_median":[],"S2_median":[],
                            "S1_dwell":[], "S2_dwell":[]} # Get sequence here from reference genome ??

                # Parse S1 and S2 reads data and add to mean and dwell time per position
                for lab, fp in (("S1", s1_fp), ("S2", s2_fp)):
                    for read in ref_dict[lab]:

                        # Move to read, save read data chunk and reset file pointer
                        fp.seek (read.byte_offset)
                        line_list = fp.read (read.byte_len).split("\n")
                        fp.seek (0)

                        # Check read_id ref_id concordance between index and data file (#6f885af6-5844-476e-9e51-133dc5617dfd	YHR055C)
                        header = line_list[0][1:].split("\t")
                        if not header[0] == read.read_id or not header[1] == read.ref_id:
                            raise NanocomporeError ("Index and data files are not matching")

                        # Extract col names from second line
                        col_names = line_list[1].split("\t")
                        line_tuple = namedtuple ("line_tuple", col_names)

                        # Parse data files kmers per kmers
                        for line in line_list[2:]:
                            lt = line_tuple (*line.split("\t"))

                            # Check if positions are in the ones found in the whitelist intervals
                            ref_pos = int(lt.ref_pos)
                            if ref_pos in ref_pos_dict:
                                # Append mean value and dwell time per position
                                ref_pos_dict[ref_pos][lab+"_median"].append (float(lt.median))
                                ref_pos_dict[ref_pos][lab+"_dwell"].append (int(lt.n_signals))

                # Filter low coverage positions
                low_cov_pos = []
                for pos, pos_dict in ref_pos_dict.items():
                    if len(pos_dict["S1_median"]) < self.__min_coverage or len(pos_dict["S2_median"]) < self.__min_coverage:
                        low_cov_pos.append (pos)
                for pos in low_cov_pos:
                    del ref_pos_dict [pos]

                if ref_pos_dict:

                    if self.__comparison_method in ["mann_whitney", "MW", "kolmogorov_smirnov", "KS","t_test", "TT"]:
                        ref_pos_dict = paired_test (
                            ref_pos_dict=ref_pos_dict,
                            method=self.__comparison_method,
                            sequence_context=self.__sequence_context,
                            min_coverage=self.__min_coverage)

                    elif self.__comparison_method == "kmean":
                        pass

                    # Add the current read details to queue
                    out_q.put ((ref_id, ref_pos_dict))
        # Add poison pill in queues
        out_q.put (None)


    def __write_output (self, out_q):
        # Get results out of the out queue and write in shelve
        with shelve.open (self.__output_db_fn, flag='n') as db:
            # Iterate over the counter queue and process items until all poison pills are found
            pbar = tqdm (total = len(self.__whitelist), unit=" Processed References", disable=self.__logLevel=="warning")
            for _ in range (self.__nthreads):
                for ref_id, ref_pos_dict in iter (out_q.get, None):
                    # Write results in a shelve db to get around multithreaded isolation
                    db [ref_id] = ref_pos_dict
                    pbar.update ()
            pbar.close()


#~~~~~~~~~~~~~~PROPERTY HELPER AND MAGIC METHODS~~~~~~~~~~~~~~#

    def __len__ (self):
        try:
            return self.__len
        except:
            with shelve.open (self.__output_db_fn, flag='r') as db:
                self.__len = len(db)
                return self.__len

    def __iter__ (self):
        with shelve.open (self.__output_db_fn, flag = "r") as db:
            for k, v in db.items():
                yield (k,v)

    def __getitem__(self, items):
        with shelve.open (self.__output_db_fn, flag = "r") as db:
            if items in db:
                return db[items]
            else:
                return None

    @property
    def ref_id_list (self):
        try:
            return self.__ref_id_list
        except:
            with shelve.open (self.__output_db_fn, flag='r') as db:
                self.__ref_id_list = list (db.keys())
                return self.__ref_id_list
