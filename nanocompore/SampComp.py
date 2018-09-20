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
from nanocompore.helper_lib import mkdir, access_file
from nanocompore.Whitelist import Whitelist
from nanocompore.NanocomporeError import NanocomporeError

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
        whitelist,
        output_db_fn,
        padj_threshold = 0.1,
        comparison_method = "kmean",
        sequence_context=0,
        nthreads = 4,
        logLevel = "info",):

        """
        s1_fn: Path to sample 1 eventalign_collapse data file
        s2_fn: Path to sample 2 eventalign_collapse data file
        whitelist: Whitelist object previously generated with nanocompore Whitelist
        output_db_fn: Path where to write the result database
        padj_threshold: Adjusted p-value threshold for reporting sites.
        comparison_method: Statistical method to compare the 2 samples signal (default kmean)
        sequence_context: Extend statistical analysis to contigous adjacent base is available
        nthreads: Number of threads (two are used for reading and writing, all the others for processing in parallel).
        logLevel: Set the log level. Valid values: warning, info, debug
        """
        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise and checks options")

        # Check args
        if not isinstance (whitelist, Whitelist):
            raise NanocomporeError("Whitelist is not valid")
        for fn in (s1_fn, s2_fn):
            if not access_file (fn):
                raise NanocomporeError("Cannot access file {}".format(fn))
        if nthreads < 3:
            raise NanocomporeError("Number of threads not valid")

        # Save private args
        self.__s1_fn = s1_fn
        self.__s2_fn = s2_fn
        self.__whitelist = whitelist
        self.__output_db_fn = output_db_fn
        self.__padj_threshold = padj_threshold
        self.__comparison_method = comparison_method
        self.__sequence_context = sequence_context
        self.__nthreads = nthreads - 2 # Remove reader and writer threads
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
        # Consumme ref_id and position_dict until empty and perform statiscical analysis

        with open (self.__s1_fn) as s1_fp, open (self.__s2_fn) as s2_fp: # More efficient to open only once the files
            for ref_id in iter (in_q.get, None):

                ref_dict = self.__whitelist[ref_id]
                position_dict = OrderedDict ()

                for interval_start, interval_end in ref_dict["interval_list"]:
                    for i in range (interval_start, interval_end+1):
                        position_dict[i] = {"S1":[], "S2":[]}

                # Parse S1 and S2 reads data and add to mean and dwell time per position
                for lab, fp in (("S1", s1_fp), ("S2", s2_fp)):
                    for read in ref_dict[lab]:

                        # Move to read save read data chunk and reset file pointer
                        fp.seek (read.byte_offset)
                        line_list = fp.read (read.byte_len).split("\n")
                        fp.seek (0)

                        # Extract info from line
                        header = line_list[1].split("\t")
                        line_tuple = namedtuple("line_tuple", header)
                        for line in line_list[2:]:
                            lt = line_tuple(*line.split("\t"))

                            # Check if positions are in the ones found in the whitelist intervals
                            ref_pos = int(lt.ref_pos)
                            if ref_pos in position_dict:
                                # Append mean value and dwell time per position
                                position_dict[ref_pos][lab].append ((float(lt.mean), int(lt.n_signals)))

                # Do stats with position_dicts
                ####### ## Add p-value per position to the position_dict #######
                ####### position_dict = tx_compare (self.__padj_threshold, self.__comparison_method, self.__sequence_context) #######
                # Add the current read details to queue
                out_q.put ((ref_id, position_dict))
        # Add poison pill in queues
        out_q.put (None)

    def __write_output (self, out_q): ############################################# Or pickle dict or flat file ...
        # Get results out of the out queue and write in shelve
        with shelve.open (self.__output_db_fn, flag='n') as db:
            # Iterate over the counter queue and process items until all poison pills are found
            pbar = tqdm (total = len(self.__whitelist), unit=" Processed References", disable=self.__logLevel=="warning")
            for _ in range (self.__nthreads):
                for ref_id, stat_dict in iter (out_q.get, None):
                    # Write results in a shelve db to get around multithreaded isolation
                    db [ref_id] = stat_dict
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
