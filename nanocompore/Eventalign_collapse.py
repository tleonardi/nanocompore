# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import multiprocessing as mp
from time import time
from collections import *
import traceback
import os

# Third party imports
from loguru import logger
import numpy as np
from tqdm import tqdm

# Local imports
from nanocompore.common import *
from nanocompore.SuperParser import SuperParser

# Disable multithreading for MKL and openBlas
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

log_level_dict = {"debug":"DEBUG", "info":"INFO", "warning":"WARNING"}
logger.remove()

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Eventalign_collapse ():

    def __init__ (self,
        eventalign_fn:str,
        outpath:str="./",
        outprefix:str="out",
        overwrite:bool = False,
        n_lines:int=None,
        nthreads:int = 3,
        progress:bool = False):
        """
        Collapse the nanopolish eventalign events at kmer level
        * eventalign_fn
            Path to a nanopolish eventalign tsv output file, or a list of file, or a regex (can be gzipped)
        * outpath
            Path to the output folder (will be created if it does exist yet)
        * outprefix
            text outprefix for all the files generated
        * overwrite
            If the output directory already exists, the standard behaviour is to raise an error to prevent overwriting existing data
            This option ignore the error and overwrite data if they have the same outpath and outprefix.
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
        self.__outpath = outpath
        self.__outprefix = outprefix
        self.__eventalign_fn = eventalign_fn
        self.__n_lines = n_lines
        self.__nthreads = nthreads - 2 # subtract 1 for reading and 1 for writing
        self.__progress = progress

        # Input file field selection typing and renaming
        self.__select_colnames = ["contig", "read_name", "position", "reference_kmer", "model_kmer", "event_length", "samples"]
        self.__change_colnames = {"contig":"ref_id" ,"position":"ref_pos", "read_name":"read_id", "samples":"sample_list", "event_length":"dwell_time"}
        self.__cast_colnames = {"ref_pos":int, "dwell_time":np.float32, "sample_list":lambda x: [float(i) for i in x.split(",")]}

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
        ps_list.append (mp.Process (target=self.__split_reads, args=(in_q, error_q)))
        for i in range (self.__nthreads):
            ps_list.append (mp.Process (target=self.__process_read, args=(in_q, out_q, error_q)))
        ps_list.append (mp.Process (target=self.__write_output, args=(out_q, error_q)))

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

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __split_reads (self, in_q, error_q):
        """
        Mono-threaded reader
        """
        logger.debug("Start reading input file(s)")

        n_reads = n_events = 0

        try:
            # Open input file with superParser
            with SuperParser(
                fn = self.__eventalign_fn,
                select_colnames = self.__select_colnames,
                cast_colnames = self.__cast_colnames,
                change_colnames = self.__change_colnames,
                n_lines=self.__n_lines) as sp:

                for l in sp:
                    n_events+=1

                    # First event exception
                    if n_events==1:
                        cur_ref_id = l["ref_id"]
                        cur_read_id = l["read_id"]
                        event_l = [l]

                    # Same read/ref group = just append to current event group
                    elif l["ref_id"] == cur_ref_id and l["read_id"] == cur_read_id:
                        event_l.append(l)

                    # If new read/ref group detected = enqueue previous event group and start new one
                    else:
                        n_reads+=1
                        in_q.put(event_l)

                        cur_ref_id = l["ref_id"]
                        cur_read_id = l["read_id"]
                        event_l = [l]

                # Last event line exception
                in_q.put(event_l)
                n_reads+=1

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.debug("Error in Reader")
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pills
        finally:
            for i in range (self.__nthreads):
                in_q.put(None)
            logger.debug("Parsed Reads:{} Events:{}".format(n_reads, n_events))

    def __process_read (self, in_q, out_q, error_q):
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
                read = Read (read_id=events_l[0]["read_id"], ref_id=events_l[0]["ref_id"])
                for event_d in events_l:
                    read.add_event(event_d)
                    n_events+=1

                # If at least one valid event found collect results at read and kmer level
                if read.n_events > 1:
                    read_res_d = read.get_read_results()
                    kmer_res_l = read.get_kmer_results()
                    out_q.put((read_res_d, kmer_res_l))
                    n_reads+=1
                    n_kmers+= len(kmer_res_l)
                    n_signals+= read.n_signals

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.error("Error in Worker. Event:{}".format(event_d))
            error_q.put (NanocomporeError(traceback.format_exc()))

        # Deal poison pill
        finally:
            logger.debug("Processed Reads:{} Kmers:{} Events:{} Signals:{}".format(n_reads, n_kmers, n_events, n_signals))
            out_q.put(None)

    def __write_output (self, out_q, error_q):
        """
        Mono-threaded Writer
        """
        logger.debug("Start writing output files")

        byte_offset = n_reads = n_kmers = 0

        # Init variables for index files
        idx_fn = os.path.join(self.__outpath, self.__outprefix+"_eventalign_collapse.tsv.idx")
        data_fn = os.path.join(self.__outpath, self.__outprefix+"_eventalign_collapse.tsv")

        try:
            # Open output files and tqdm progress bar
            with open (data_fn, "w") as data_fp, open (idx_fn, "w") as idx_fp, tqdm (unit=" reads", disable=not self.__progress) as pbar:

                # Iterate over out queue until nthread poison pills are found
                for _ in range (self.__nthreads):
                    for read_res_d, kmer_res_l in iter (out_q.get, None):
                        n_reads+=1

                        # Define file header from first read and first kmer
                        if byte_offset is 0:
                            idx_header_list = list(read_res_d.keys())+["byte_offset","byte_len"]
                            idx_header_str = "\t".join(idx_header_list)
                            data_header_list = list(kmer_res_l[0].keys())
                            data_header_str = "\t".join(data_header_list)

                            # Write index file header
                            idx_fp.write ("{}\n".format(idx_header_str))

                        # Write data file header
                        byte_len = 0
                        header_str = "#{}\t{}\n{}\n".format(read_res_d["read_id"], read_res_d["ref_id"], data_header_str)
                        data_fp.write(header_str)
                        byte_len+=len(header_str)

                        # Write kmer data matching data field order
                        for kmer in kmer_res_l:
                            n_kmers+=1
                            data_str = "\t".join([str(kmer[f]) for f in data_header_list])+"\n"
                            data_fp.write(data_str)
                            byte_len+=len(data_str)

                        # Add byte
                        read_res_d["byte_offset"] = byte_offset
                        read_res_d["byte_len"] = byte_len-1
                        idx_str = "\t".join([str(read_res_d[f]) for f in idx_header_list])
                        idx_fp.write("{}\n".format(idx_str))

                        # Update pbar
                        byte_offset+=byte_len
                        pbar.update(1)

                # Flag last line
                data_fp.write ("#\n")

        # Manage exceptions and add error trackback to error queue
        except Exception:
            logger.error("Error in Writer")
            error_q.put (NanocomporeError(traceback.format_exc()))

        finally:
            logger.debug("Written Reads:{} Kmers:{}".format(n_reads, n_kmers))
            logger.info ("Output reads written:{}".format(n_reads))
            # Kill error queue with poison pill
            error_q.put(None)

#~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASSES~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Read ():
    """Helper class representing a single read"""

    def __init__ (self, read_id, ref_id):
        self.read_id = read_id
        self.ref_id = ref_id
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
        d["kmers"] = self.ref_end-self.ref_start+1
        d["missing_kmers"] = d["kmers"] - len(self.kmer_l)
        d["NNNNN_kmers"]=0
        d["mismatch_kmers"]=0
        d["valid_kmers"]=0
        for k in self.kmer_l:
            d[k.status+"_kmers"]+=1
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

class Kmer ():
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
        sample_array = np.array(self.sample_list, dtype=np.float32)
        d = OrderedDict()
        d["ref_pos"] = self.ref_pos
        d["ref_kmer"] = self.ref_kmer
        d["num_events"] = self.num_events
        d["num_signals"] = self.num_signals
        d["dwell_time"] = self.dwell_time
        d["NNNNN_dwell_time"] = self.NNNNN_dwell_time
        d["mismatch_dwell_time"] = self.mismatch_dwell_time
        d["status"] = self.status
        d["median"] = np.median(sample_array)
        d["mad"] = np.median(np.abs(sample_array-d["median"]))

        return d
