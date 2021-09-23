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

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Eventalign_collapse ():

    def __init__(self,
                 sample_dict:dict,
                 output_dir:str,
                 master_db:str = "nanocompore.db",
                 output_subdirs:int = 100,
                 overwrite:bool = False,
                 n_lines:int = None,
                 nthreads:int = 3):
        """
        Collapse the nanopolish eventalign events at kmer level
        * sample_dict
            Dictionary with input file and sample information:
            {condition1: [(file1, sample1), (file2, sample2), ...], condition2: [...]}
        * output_dir
            Path to the output directory
        * output_subdirs
            Distribute output files into this many subdirectories
        * overwrite
            Overwrite an existing output file?
        * n_lines
            Maximum number of read to parse.
        * nthreads
            Number of threads (need at least one each for the main process, reading, and writing).
        """
        logger.info("Checking and initialising Eventalign_collapse")

        # Save init options in dict for later
        log_init_state(loc=locals())

        # Check threads number
        if nthreads < 3:
            raise NanocomporeError("Minimal required number of threads is 3")

        # Check sample information
        self._sample_dict = sample_dict
        self.__check_sample_dict()

        # Save args to self values
        self._output_dir = output_dir
        self._master_db_path = os.path.join(output_dir, master_db)
        self._output_subdirs = max(1, output_subdirs)
        self._overwrite = overwrite
        self._n_lines = n_lines
        self._nthreads = nthreads - 1 # subtract 1 for main process

        # Input file field selection typing and renaming
        self._select_colnames = ["contig", "read_name", "position", "reference_kmer", "model_kmer", "event_length", "samples"]
        self._change_colnames = {"contig": "ref_id", "position": "ref_pos", "read_name": "read_id", "samples": "sample_list", "event_length": "dwell_time"}
        self._cast_colnames = {"ref_pos": int, "dwell_time": np.float32, "sample_list": lambda x: [float(i) for i in x.split(",")]}


    def __check_sample_dict(self):
        # expected format:
        # {condition1: [(file1, sample1), (file2, sample2), ...], condition2: [...]}
        try:
            conditions = list(self._sample_dict.keys())
            if len(conditions != 2):
                raise NanocomporeError(f"Expected two experimental conditions, found {len(conditions)}: {', '.join(conditions)}")
            samples = set([])
            paths = set([])
            for condition, sample_list in self._sample_dict.items():
                if len(sample_list) == 0:
                    raise NanocomporeError(f"No input file/sample information for condition '{condition}'")
                for path, sample in sample_list:
                    if sample in samples:
                        raise NanocomporeError(f"Sample labels must be unique. Found duplicate label: {sample}.")
                    samples.add(sample)
                    if not os.path.isfile(path):
                        raise NanocomporeError(f"Input file not found: {path}")
                    if path in paths:
                        raise NanocomporeError(f"Found duplicate input file: {path}.")
                    paths.add(path)
        except:
            logger.error("Input file/sample information failed to validate")
            raise


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

        in_q = mp.Queue() # queue for input files
        out_qs = [] # queues for reads (stratified by transcript accession)
        error_q = mp.Queue() # queue for errors
        master_db_lock = mp.Lock() # lock for write access to master DB

        logger.info("Creating master database, storing sample information and parameters")
        db_create_mode = DBCreateMode.OVERWRITE if self._overwrite else DBCreateMode.CREATE_MAYBE
        n_samples = 0
        with DataStore_master(self._master_db_path, db_create_mode) as db:
            for condition, sample_list in self._sample_dict.items():
                for path, sample in sample_list:
                    sample_id = db.store_sample(sample, path, condition)
                    in_q.put((sample_id, path))
                    n_samples += 1
            db.store_parameters("EC", output_subdirs=self._output_subdirs, n_lines=self._n_lines)

        # Define processes
        ps_list = []
        # Do we have enough threads to allocate one to each input file?
        if self._nthreads < n_samples * 2: # no - assign half the threads to reading/writing
            self._n_readers = self._nthreads // 2
        else: # yes - use one reader thread per input, the rest writers
            self._n_readers = n_samples
        for i in range(self._n_readers):
            in_q.put(None) # one "poison pill" to be consumed by each reader thread
            ps_list.append(mp.Process(target=self.__split_input, args=(in_q, out_qs, error_q, master_db_lock)))
        self._n_writers = self._nthreads - self._n_readers
        logger.info(f"Using {n_readers} reader and {n_writers} writer processes")
        # Each writer thread is associated with a specific queue (to avoid access conflicts on transcript DBs):
        for i in range(self._n_writers):
            out_qs.append(mp.Queue(maxsize=100))
            ps_list.append(mp.Process(target=self.__process_reads, args=(out_qs[i], error_q)))
        q_list = [in_q, error_q] + out_qs

        logger.info("Starting data processing")
        try:
            # Start all processes
            for ps in ps_list:
                ps.start()
            # Monitor error queue
            for _ in range(self._n_writers):
                for trace in iter(error_q.get, None):
                    logger.error("Error caught from error_q")
                    raise NanocomporeError(trace)
            # Soft processes and queues stopping
            for ps in ps_list:
                ps.join()
            for q in q_list:
                q.close()

        # Catch error, kill all processed and reraise error
        except Exception as E:
            logger.error("An error occured. Killing all processes and closing queues")
            try:
                for ps in ps_list:
                    ps.terminate()
                for q in q_list:
                    q.close()
            except:
                logger.error("An error occured while trying to kill processes")
            raise E


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __handle_transcript(self, accession, tx_info, n_qs, master_db_lock):
        """Deal with a transcript accession found by '__split_input'"""
        if accession not in tx_info:
            # Calculate two bin numbers derived from the hash:
            q_index, subdir = get_hash_bin(accession, (n_qs, self._output_subdirs))
            # Add transcript to master DB (or get ID if it already exists):
            with master_db_lock, DataStore_master(self._master_db_path) as db:
                tx_id = db.store_transcript(accession, str(subdir))
            logger.debug(f"Added new transcript to DB: {accession}")
            tx_info[accession] = (tx_id, q_index, subdir)
        return tx_info[accession]


    def __split_input(self, in_q, out_qs, error_q, master_db_lock):
        """Split input data into reads and enqueue them for processing/output"""
        tx_info = {} # information on transcripts (keyed by accession)
        try:
            for sample_id, path in iter(in_q.get, None):
                n_reads = n_events = 0
                logger.info(f"Reading input file: {path}")
                with SuperParser(path,
                                 select_colnames=self._select_colnames,
                                 cast_colnames=self._cast_colnames,
                                 change_colnames=self._change_colnames,
                                 n_lines=self._n_lines) as sp:
                    # First line/event - initialise
                    line = next(iter(sp))
                    # TODO: read ID should be unique, so no need to check transcript - correct?
                    cur_read_id = line["read_id"]
                    events = [line]
                    n_events = 1
                    # All following lines
                    for line in sp:
                        n_events += 1
                        # Same read - just append to current event group
                        if line["read_id"] == cur_read_id:
                            events.append(line)
                        # New read/ref group detected - enqueue previous event group and start new one
                        else:
                            n_reads += 1
                            tx_id, q_index, subdir = self.__handle_transcript(line["ref_id"], tx_info, len(out_qs),
                                                                              master_db_lock)
                            out_qs[q_index].put((sample_id, tx_id, subdir, events))
                            cur_read_id = line["read_id"]
                            events = [line]
                    # Last event/line
                    n_reads += 1
                    tx_id, q_index, subdir = self.__handle_transcript(line["ref_id"], tx_info, len(out_qs),
                                                                      master_db_lock)
                    out_qs[q_index].put((sample_id, tx_id, subdir, events))
                logger.debug(f"Parsed {n_reads} reads, {n_events} events from {len(tx_info)} transcripts in input file '{path}'")
        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in reader process")
            error_q.put(NanocomporeError(traceback.format_exc()))
        # Deal poison pills
        finally:
            for out_q in out_qs:
                out_q.put(None)


    def __process_reads(self, out_q, error_q):
        """Process reads from output queue and write to transcript database"""
        logger.debug("Starting processing reads")
        try:
            n_reads = n_kmers = n_events = n_signals = 0
            for _ in range(self._n_readers): # every reader is going to put items on the queue
                # Take one event list corresponding to one read from the list
                for sample_id, tx_id, subdir, events in iter(out_q.get, None):
                    # Create an empty Read object and fill it with event lines
                    # events aggregation at kmer level is managed within the object
                    read = Read(read_id=events[0]["read_id"], ref_id=events[0]["ref_id"], sample_name=sample_id)
                    for event in events:
                        read.add_event(event)
                        n_events += 1
                    # If at least one valid event found collect results at read and kmer level
                    if read.n_events > 1:
                        n_reads += 1
                        n_kmers += len(read.kmer_l)
                        n_signals += read.n_signals
                        # Write read to transcript database:
                        db_path = os.path.join(self._output_dir, str(subdir), read.ref_id + ".db")
                        with DataStore_transcript(db_path, read.ref_id, tx_id, DBCreateMode.CREATE_MAYBE) as db:
                            db.store_read(read, sample_id)
        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.error("Error in writer process.")
            error_q.put(NanocomporeError(traceback.format_exc()))
        # Deal poison pill
        finally:
            logger.debug(f"Processed {n_reads} reads, {n_kmers} kmers, {n_events} events, {n_signals} signals")
            error_q.put(None)


#~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASSES~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Read:
    """Helper class representing a single read"""

    def __init__(self, read_id, ref_id, sample_name):
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

    def add_event(self, event_d):
        self.n_events += 1
        self.n_signals += len(event_d["sample_list"])
        self.dwell_time += event_d["dwell_time"]
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
    def kmers_status(self):
        d = OrderedDict()
        d["kmers"] = self.ref_end - self.ref_start + 1
        d["missing_kmers"] = d["kmers"] - len(self.kmer_l)
        d["NNNNN_kmers"] = 0
        d["mismatch_kmers"] = 0
        d["valid_kmers"] = 0
        for k in self.kmer_l:
            d[k.status + "_kmers"] += 1
        return d

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
