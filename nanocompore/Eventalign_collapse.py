# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import multiprocessing as mp
from time import time
from collections import *
import traceback
import os
import cProfile as profile
import statistics

# Third party imports
from loguru import logger
import numpy as np
from tqdm import tqdm

# Local imports
from nanocompore.common import *
from nanocompore.SuperParser import SuperParser
from nanocompore.DataStore import DataStore_master, DataStore_transcript, DBCreateMode

# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

log_level_dict = {"debug": "DEBUG", "info": "INFO", "warning": "WARNING"}
#logger.remove()

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Eventalign_collapse ():

    def __init__(self,
                 eventalign_path:str,
                 sample_name:str,
                 condition:str,
                 output_dir:str,
                 output_subdirs:int = 100,
                 overwrite:bool = False,
                 n_lines:int = None,
                 nthreads:int = 3,
                 progress:bool = False):
        # TODO: is 'overwrite' a useful option, as data from multiple samples needs to be accumulated in the same DB?
        """
        Collapse the nanopolish eventalign events at kmer level
        * eventalign_path
            Path to a nanopolish eventalign tsv output file, or a list of files, or a regex (can be gzipped)
        * sample_name
            The name of the sample being processed
        * condition
            Condition label of the sample
        * output_dir
            Path to the output directory
        * output_subdirs
            Distribute output files into this many subdirectories
        * overwrite
            Overwrite an existing output file?
        * n_lines
            Maximum number of read to parse.
        * nthreads
            Number of threads (two are used for reading and writing, all the others for parallel processing).
        * progress
            Display a progress bar during execution
        """
        logger.info("Checking and initialising Eventalign_collapse")

        # Save init options in dict for later
        log_init_state(loc=locals())

        # Check threads number
        if nthreads < 3:
            raise NanocomporeError("Minimal required number of threads >=3")

        # Save args to self values
        self._eventalign_path = eventalign_path
        self._sample_name = sample_name
        self._condition = condition
        self._output_dir = output_dir
        self._output_subdirs = max(1, output_subdirs)
        self._overwrite = overwrite
        self._n_lines = n_lines
        self._nthreads = nthreads - 2 # subtract 1 for reading and 1 for writing
        self._progress = progress

        # Input file field selection typing and renaming
        self._select_colnames = ["contig", "read_name", "position", "reference_kmer", "model_kmer", "event_length", "samples"]
        self._change_colnames = {"contig": "ref_id", "position": "ref_pos", "read_name": "read_id", "samples": "sample_list", "event_length": "dwell_time"}
        self._cast_colnames = {"ref_pos":int, "dwell_time":np.float32, "sample_list":lambda x: [float(i) for i in x.split(",")]}


    def __call__(self):
        """
        Run the analysis
        """
        logger.info("Creating output subdirectories")
        for i in range(self._output_subdirs):
            subdir = os.path.normpath(os.path.join(self._output_dir, str(i)))
            if self._overwrite and os.path.exists(subdir):
                # TODO: what if 'subdir' is a file, not a directory?
                shutil.rmtree(subdir)
            os.makedirs(subdir, exist_ok=True)

        logger.info("Starting data processing")
        # Init Multiprocessing variables
        in_q = mp.Queue(maxsize = 100)
        out_q = mp.Queue(maxsize = 100)
        error_q = mp.Queue()

        # Define processes
        ps_list = []
        ps_list.append(mp.Process(target=self.__split_reads, args=(in_q, error_q)))
        for i in range (self._nthreads):
            ps_list.append(mp.Process(target=self.__process_read, args=(in_q, out_q, error_q)))
        # TODO: Check that sample_name does not exist already in DB
        ps_list.append(mp.Process(target=self.__write_output, args=(out_q, error_q)))

        # Start processes and monitor error queue
        try:
            # Start all processes
            for ps in ps_list:
                ps.start()

            # Monitor error queue
            for tb in iter(error_q.get, None):
                logger.error("Error caught from error_q")
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


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __split_reads(self, in_q, error_q):
        """
        Mono-threaded reader
        """
        logger.debug("Start reading input file(s)")

        n_reads = n_events = 0

        try:
            # Open input file with superParser
            # TODO: benchmark performance compared to csv.DictReader (std. lib.)
            with SuperParser(
                fn = self._eventalign_path,
                select_colnames = self._select_colnames,
                cast_colnames = self._cast_colnames,
                change_colnames = self._change_colnames,
                n_lines = self._n_lines) as sp:

                # First line/event - initialise
                l = next(iter(sp))
                # TODO: read ID should be unique, so no need to check transcript - correct?
                cur_read_id = l["read_id"]
                event_l = [l]
                n_events = 1

                # All following lines
                for l in sp:
                    n_events += 1
                    # Same read = just append to current event group
                    if l["read_id"] == cur_read_id:
                        event_l.append(l)
                    # If new read/ref group detected = enqueue previous event group and start new one
                    else:
                        n_reads += 1
                        in_q.put(event_l)
                        cur_read_id = l["read_id"]
                        event_l = [l]

                # Last event/line
                in_q.put(event_l)
                n_reads += 1

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in Reader")
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pills
        finally:
            for i in range (self._nthreads):
                in_q.put(None)
            logger.debug("Parsed Reads:{} Events:{}".format(n_reads, n_events))


    def __process_read(self, in_q, out_q, error_q):
        """
        Multi-threaded workers collapsing events at kmer level
        """
        logger.debug("Starting processing reads")
        try:
            n_reads = n_kmers = n_events = n_signals = 0

            # Take on event list corresponding to one read from the list
            for events_l in iter(in_q.get, None):

                # Create an empty Read object and fill it with event lines
                # events aggregation at kmer level is managed withon the object
                read = Read(read_id=events_l[0]["read_id"], ref_id=events_l[0]["ref_id"],
                            sample_name=self._sample_name)
                for event_d in events_l:
                    read.add_event(event_d)
                    n_events += 1

                # If at least one valid event found collect results at read and kmer level
                if read.n_events > 1:
                    read_res_d = read.get_read_results()
                    kmer_res_l = read.get_kmer_results()
                    out_q.put(read)
                    n_reads += 1
                    n_kmers += len(kmer_res_l)
                    n_signals += read.n_signals

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.error("Error in Worker.")
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pill
        finally:
            logger.debug("Processed Reads:{} Kmers:{} Events:{} Signals:{}".format(n_reads, n_kmers, n_events, n_signals))
            out_q.put(None)


    def __write_output(self, out_q, error_q):
        """
        Single-threaded writer
        """
        logger.debug("Start writing output to DB")

        # pr = profile.Profile()
        # pr.enable()
        n_reads = 0
        db_create_mode = DBCreateMode.OVERWRITE if self._overwrite else DBCreateMode.CREATE_MAYBE
        master_db_path = os.path.join(self._output_dir, "eventalign_collapse.db")
        try:
            with DataStore_master(master_db_path, db_create_mode) as master, tqdm(unit=" reads") as pbar:
                sample_id = master.store_sample(self._sample_name, self._condition)
                # Iterate over out queue until nthread poison pills are found
                for _ in range(self._nthreads):
                    for read in iter(out_q.get, None):
                        n_reads += 1
                        # Which database file to use? - depends on transcript:
                        accession = read.ref_id
                        # To distribute transcripts evenly over subdirectories, use (part of) the hash value:
                        # Don't use built-in 'hash' function!
                        # Values are salted using a random seed, so not reproducible between Python runs!
                        subdir = str(get_hash_bin(accession, self._output_subdirs))
                        tx_id = master.store_transcript(accession, subdir)
                        # TODO: what if accession contains characters that aren't allowed in paths?
                        db_path = os.path.join(self._output_dir, subdir, accession + ".db")
                        with DataStore_transcript(db_path, accession, tx_id, DBCreateMode.CREATE_MAYBE) as db:
                            db.store_read(read, sample_id)
                        pbar.update(1)
        except Exception:
            logger.error("Error adding read to DB")
            error_q.put(NanocomporeError(traceback.format_exc()))

        finally:
            logger.info ("Output reads written:{}".format(n_reads))
            # Kill error queue with poison pill
            error_q.put(None)
            # pr.disable()
            # pr.dump_stats("prof")


#~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASSES~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Read:
    """Helper class representing a single read"""

    def __init__ (self, read_id, ref_id, sample_name):
        self.read_id = read_id
        self.ref_id = ref_id
        self.sample_name = sample_name
        self.ref_start = None
        self.ref_end = None
        self.kmer_l = [Kmer()]
        self.n_events = 0
        self.n_signals = 0
        self.dwell_time = 0

    def __repr__(self):
        s = "Read instance\n"
        s+="\tread_id: {}\n".format(self.read_id)
        s+="\tref_id: {}\n".format(self.ref_id)
        s+="\tref_start: {}\n".format(self.ref_start)
        s+="\tref_end: {}\n".format(self.ref_end)
        s+="\tn_events: {}\n".format(self.n_events)
        s+="\tdwell_time: {}\n".format(self.dwell_time)
        for status, count in self.kmers_status.items():
            s+="\t{}: {}\n".format(status, count)
        return s

    def add_event (self, event_d):
        self.n_events+=1
        self.n_signals+=len(event_d["sample_list"])
        self.dwell_time+=event_d["dwell_time"]

        # First event
        if not self.ref_end:
            self.ref_start = event_d["ref_pos"]
            self.ref_end = event_d["ref_pos"]

        # Position offset = Move to next position
        if event_d["ref_pos"] > self.ref_end:
            self.kmer_l.append(Kmer())
            self.ref_end = event_d["ref_pos"]
            self.ref_end = event_d["ref_pos"]

        # Update current kmer
        self.kmer_l[-1].add_event(event_d)

    @property
    def kmers_status (self):
        d = OrderedDict()
        d["kmers"] = self.ref_end - self.ref_start + 1
        d["missing_kmers"] = d["kmers"] - len(self.kmer_l)
        d["NNNNN_kmers"] = 0
        d["mismatch_kmers"] = 0
        d["valid_kmers"] = 0
        for k in self.kmer_l:
            d[k.status + "_kmers"] += 1
        return d

    def get_read_results (self):
        d = OrderedDict()
        d["ref_id"] = self.ref_id
        d["ref_start"] = self.ref_start
        d["ref_end"] = self.ref_end
        d["read_id"] = self.read_id
        d["num_events"] = self.n_events
        d["num_signals"] = self.n_signals
        d["dwell_time"] = self.dwell_time
        for status, count in self.kmers_status.items():
            d[status] = count
        return d

    def get_kmer_results (self):
        l = [kmer.get_results() for kmer in self.kmer_l]
        return l

class Kmer:
    """Helper class representing a single kmer"""

    def __init__ (self):
        self.ref_pos = None
        self.ref_kmer = None
        self.num_events = 0
        self.num_signals = 0
        self.dwell_time = 0.0
        self.NNNNN_dwell_time = 0.0
        self.mismatch_dwell_time = 0.0
        self.sample_list = []

    def __repr__(self):
        s = "Kmer instance\n"
        s+="\tref_pos: {}\n".format(self.ref_pos)
        s+="\tref_kmer: {}\n".format(self.ref_kmer)
        s+="\tnum_events: {}\n".format(self.num_events)
        s+="\tdwell_time: {}\n".format(self.dwell_time)
        s+="\tNNNNN_dwell_time: {}\n".format(self.NNNNN_dwell_time)
        s+="\tmismatch_dwell_time: {}\n".format(self.mismatch_dwell_time)
        s+="\tsample_list: {}\n".format(" ".join([str(i) for i in self.sample_list]))
        return s

    def add_event (self, event_d):
        # First event
        if not self.ref_pos:
            self.ref_pos = event_d["ref_pos"]
            self.ref_kmer = event_d["reference_kmer"]

        # Update kmer values
        self.num_events+=1
        self.num_signals+=len(event_d["sample_list"])
        self.dwell_time+=event_d["dwell_time"]
        if event_d["model_kmer"] == "NNNNN":
            self.NNNNN_dwell_time+=event_d["dwell_time"]
        elif event_d["model_kmer"] != self.ref_kmer:
            self.mismatch_dwell_time+=event_d["dwell_time"]
        self.sample_list.extend(event_d["sample_list"])

    @property
    def status (self):
        """Define kmer status depending on mismatching or NNNNN events"""
        # If no NNNNN nor mismatch events then the kmer the purely valid
        if not self.NNNNN_dwell_time and not self.mismatch_dwell_time:
            return "valid"
        # Otherwise return status corresponding to longuest invalid dwell time
        elif self.NNNNN_dwell_time > self.mismatch_dwell_time:
            return "NNNNN"
        else:
            return "mismatch"

    def get_results(self):
        d = OrderedDict()
        d["ref_pos"] = self.ref_pos
        d["ref_kmer"] = self.ref_kmer
        d["num_events"] = self.num_events
        d["num_signals"] = self.num_signals
        d["dwell_time"] = self.dwell_time
        d["NNNNN_dwell_time"] = self.NNNNN_dwell_time
        d["mismatch_dwell_time"] = self.mismatch_dwell_time
        d["status"] = self.status
        d["median"] = statistics.median(self.sample_list)
        d["mad"] = statistics.median([abs(i - d["median"]) for i in self.sample_list])
        return d
