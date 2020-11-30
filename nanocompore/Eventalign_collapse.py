# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import multiprocessing as mp
from time import time
from collections import *
import traceback
import datetime
import os

# Third party imports
from loguru import logger
import numpy as np
from tqdm import tqdm

# Local imports
from NanopolishComp.common import *
from NanopolishComp import __version__ as package_version
from NanopolishComp import __name__ as package_name

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
        max_reads:int=None,
        write_samples:bool=False,
        stat_fields:list=["mean", "median", "num_signals"],
        nthreads:int = 3,
        log_level:str = "info"):
        """
        Collapse the nanopolish eventalign output by kmers rather that by events.
        kmer level statistics (mean, median, std, mad) are only computed if nanopolish is run with --samples option
        * eventalign_fn
            Path to a nanopolish eventalign tsv output file.
        * outpath
            Path to the output folder (will be created if it does exist yet)
        * outprefix
            text outprefix for all the files generated
        * max_reads
            Maximum number of read to parse. 0 to deactivate (default = 0)
        * write_samples
            If given, will write the raw sample if nanopolish eventalign was ran with --samples option
        * stat_fields
            List of statistical fields to compute if nanopolish eventalign was ran with --samples option.
            Valid values = "mean", "std", "median", "mad", "num_signals"
        * threads
            Total number of threads. 1 thread is used for the reader and 1 for the writer (default = 4)
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """

        # Save init options in dict for later
        kwargs = locals()
        option_d = OrderedDict()
        option_d["package_name"] = package_name
        option_d["package_version"] = package_version
        option_d["timestamp"] = str(datetime.datetime.now())
        for i, j in kwargs.items():
            if not i in ["self"]:
                option_d[i]=j

        # Check if output folder already exists
        try:
            mkdir(fn=outpath, exist_ok=overwrite)
        except NanocomporeError:
            raise NanocomporeError("Could not create the output folder. Try using `overwrite` option or use another directory")

        # Write init options to log file
        log_fn = os.path.join(outpath, outprefix+"Eventalign_collapse.log")
        with open(log_fn, "w") as log_fp:
            json.dump(option_d, log_fp, indent=2)

        # Set logging level
        logger.add(sys.stderr, format="{time} {level} - {process.name} | {message}", enqueue=True, level=log_level_dict.get(log_level, "WARNING"))
        logger.add(log_fn, format="{time} {level} - {process.name} | {message}", enqueue=True, level="TRACE")
        logger.info("Initialising Eventalign_collapse and checking options")

        # Try to read input file if not a stream
        logger.debug("\tTesting input file readability")
        if eventalign_fn != 0 and not file_readable (eventalign_fn):
            raise IOError ("Cannot read input file")

        # Check at least 3 threads
        if nthreads < 3:
            raise NanocomporeError("The minimum number of threads is 3")

        logger.debug("\tChecking if stat_fields names are valid")
        for field in stat_fields:
            valid_fields = ["mean", "std", "median", "mad", "num_signals"]
            if not field in valid_fields:
                raise ValueError ("Invalid value in stat_field {}. Valid entries = {}".format(field, ",".join(valid_fields)))

        # Save args to self values
        self.outpath = outpath
        self.outprefix = outprefix
        self.eventalign_fn = eventalign_fn
        self.threads = threads-2 # Remove 2 threads for read and write
        self.max_reads = max_reads
        self.write_samples = write_samples
        self.stat_fields = stat_fields

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
        ps_list.append (mp.Process (target=self._split_reads, args=(in_q, error_q)))
        for i in range (self.threads):
            ps_list.append (mp.Process (target=self._process_read, args=(in_q, out_q, error_q, i+1)))
        ps_list.append (mp.Process (target=self._write_output, args=(out_q, error_q)))

        # Start processes and monitor error queue
        try:
            # Start all processes
            for ps in ps_list:
                ps.start()
            # Monitor error queue
            for tb in iter(error_q.get, None):
                logger.trace("Error caught from error_q")
                raise NanocomporeError(tb)

        # Catch error and reraise it
        except(BrokenPipeError, KeyboardInterrupt, NanocomporeError) as E:
            logger.error("An error occured. Killing all processes\n")
            raise E

        finally:
            # Soft processes stopping
            for ps in ps_list:
                ps.join()

            # Hard failsafe processes killing
            for ps in ps_list:
                if ps.exitcode == None:
                    ps.terminate()

    def __repr__ (self):
        m = "General options:\n"
        m+=dict_to_str(self.option_d)
        return m

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def _split_reads (self, in_q, error_q):
        """
        Mono-threaded reader
        """
        logger.debug("\t[split_reads] Start reading input file/stream")
        try:
            # Open input file or stdin if 0
            with open (self.eventalign_fn) as fp:

                # Get header line and extract corresponding index
                input_header = fp.readline().rstrip().split("\t")
                if not input_header:
                    raise NanopolishCompError ("Input file/stream is empty")

                idx = self._get_field_idx (input_header)
                n_reads = 0

                # First data line exception
                read_l = []
                event_l = fp.readline().rstrip().split("\t")
                cur_read_id = event_l[idx["read_id"]]
                cur_ref_id = event_l[idx["ref_id"]]
                event_d = self._event_list_to_dict (event_l, idx)
                read_l.append (event_d)

                for line in fp:
                    # Early ending if required
                    if self.max_reads and n_reads == self.max_reads:
                        break
                    # Get event line
                    event_l = line.rstrip().split("\t")
                    read_id = event_l[idx["read_id"]]
                    ref_id = event_l[idx["ref_id"]]
                    event_d = self._event_list_to_dict (event_l, idx)

                    # Line correspond to the same ids
                    if read_id != cur_read_id or ref_id != cur_ref_id:
                        in_q.put ((cur_read_id, cur_ref_id, read_l))
                        n_reads+=1
                        read_l = []
                        cur_read_id = read_id
                        cur_ref_id = ref_id

                    # In any case extend list corresponding to current read_id/ref_id
                    read_l.append(event_d)

                # Last data line exception
                in_q.put ((cur_read_id, cur_ref_id, read_l))
                n_reads+=1

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            for i in range (self.threads):
                in_q.put(None)
            logger.debug("\t[split_reads] Done")

    def _process_read (self, in_q, out_q, error_q, pid):
        """
        Multi-threaded workers
        """
        logger.debug("\t[process_read {}] Starting processing reads".format(pid))
        try:
            # Collapse event at kmer level
            for read_id, ref_id, read_l in iter(in_q.get, None):

                # Write read header to str
                read_str = "#{}\t{}\n".format(read_id, ref_id)
                read_str+= "{}\n".format (self._make_ouput_header(event_d=read_l[0]))

                # Init values for first kmer
                kmer_d = self._init_kmer_dict(event_d=read_l[0])

                # Init read dictionary
                read_d = OrderedDict ()
                read_d["read_id"] = read_id
                read_d["ref_id"] = ref_id
                read_d["dwell_time"] = 0.0
                read_d["kmers"] = 0
                read_d["NNNNN_kmers"] = 0
                read_d["mismatch_kmers"] = 0
                read_d["missing_kmers"] = 0
                read_d["ref_start"] = kmer_d["ref_pos"]

                # Iterate over the rest of the lines
                for event_d in read_l [1:]:
                    pos_offset = event_d["ref_pos"]-kmer_d["ref_pos"]

                    # Same position = update current kmer
                    if pos_offset == 0:
                        kmer_d = self._update_kmer_dict(kmer_d=kmer_d, event_d=event_d)

                    # New position = write previous kmer and start new one
                    else:
                        # Update read counter
                        read_d["dwell_time"] += kmer_d["dwell_time"]
                        if kmer_d["NNNNN_dwell_time"]:
                            read_d["NNNNN_kmers"] += 1
                        if kmer_d ["mismatch_dwell_time"]:
                            read_d["mismatch_kmers"] += 1
                        if pos_offset >=2:
                            read_d["missing_kmers"] += (pos_offset-1)
                        read_d["kmers"] += 1
                        # Converts previous kmer to str and init new kmer
                        read_str += "{}\n".format(self._kmer_dict_to_str(kmer_d=kmer_d))
                        kmer_d = self._init_kmer_dict(event_d=event_d)

                # Last read_d update
                read_d["dwell_time"] += kmer_d["dwell_time"]
                if kmer_d ["NNNNN_dwell_time"]:
                    read_d["NNNNN_kmers"] += 1
                if kmer_d ["mismatch_dwell_time"]:
                    read_d["mismatch_kmers"] += 1
                if pos_offset >=2:
                    read_d["missing_kmers"] += (pos_offset-1)
                read_d["ref_end"] = kmer_d["ref_pos"]+1
                read_d["kmers"] += 1

                # Last kmer
                read_str += "{}\n".format(self._kmer_dict_to_str(kmer_d=kmer_d))

                # Add the current read details to queue
                out_q.put((read_d, read_str))

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            logger.debug("\t[process_read {}] Done".format(pid))
            out_q.put(None)

    def _write_output (self, out_q, error_q):
        """
        Mono-threaded Writer
        """
        logger.debug("\t[write_output] Start rwriting output")

        byte_offset = n_reads = 0
        t = time()

        try:
            # Open output files
            data_fn = os.path.join(self.outpath, self.outprefix+"_eventalign_collapse.tsv")
            idx_fn = os.path.join(self.outpath, self.outprefix+"_eventalign_collapse.tsv.idx")
            with open (data_fn, "w") as data_fp,\
                 open (idx_fn, "w") as idx_fp,\
                 tqdm (unit=" reads", mininterval=0.1, smoothing=0.1, disable=logger.level>=30) as pbar:

                idx_fp.write ("ref_id\tref_start\tref_end\tread_id\tkmers\tdwell_time\tNNNNN_kmers\tmismatch_kmers\tmissing_kmers\tbyte_offset\tbyte_len\n")

                n_reads = 0
                for _ in range (self.threads):
                    for (read_d, read_str) in iter (out_q.get, None):
                        byte_len = len(read_str)

                        data_fp.write (read_str)
                        idx_fp.write ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (
                            read_d["ref_id"],
                            read_d["ref_start"],
                            read_d["ref_end"],
                            read_d["read_id"],
                            read_d["kmers"],
                            read_d["dwell_time"],
                            read_d["NNNNN_kmers"],
                            read_d["mismatch_kmers"],
                            read_d["missing_kmers"],
                            byte_offset,
                            byte_len-1))

                        byte_offset += byte_len
                        n_reads += 1
                        if logger.level<30:
                            pbar.update(1)

                # Flag last line
                data_fp.write ("#\n")

            # Open log file
            log_fn = os.path.join(self.outpath, self.outprefix+"_eventalign_collapse.log")
            with open (log_fn, "w") as log_fp:
                log_fp.write (str(self))

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            logger.debug("\t[write_output] Done")
            logger.warning ("[Eventalign_collapse] total reads: {} [{} reads/s]\n".format(n_reads, round (n_reads/(time()-t), 2)))
            error_q.put(None)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER PRIVATE METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def _get_field_idx (self, input_header):
        """"""
        # Get index of fields to fetch
        idx = OrderedDict()
        idx["ref_id"] = input_header.index ("contig")
        if "read_name" in input_header:
            idx["read_id"] = input_header.index ("read_name")
        elif "read_index" in input_header:
            idx["read_id"] = input_header.index ("read_index")
        idx["ref_pos"] = input_header.index ("position")
        idx["ref_kmer"] = input_header.index ("reference_kmer")
        idx["mod_kmer"] = input_header.index ("model_kmer")
        idx["event_len"] = input_header.index ("event_length")
        # Facultative field start and end index
        if "start_idx" in input_header and "end_idx" in input_header:
            idx["start_idx"] = input_header.index ("start_idx")
            idx["end_idx"] = input_header.index ("end_idx")
        # Facultative field samples
        if "samples" in input_header:
            idx["samples"] = input_header.index ("samples")
        return idx

    def _event_list_to_dict (self, event_l, idx):
        """Get interesting fields from event list and cast in appropriate type"""
        event_d = OrderedDict()
        event_d["ref_pos"] = int(event_l[idx["ref_pos"]])
        event_d["ref_kmer"] = event_l[idx["ref_kmer"]]
        event_d["mod_kmer"] = event_l[idx["mod_kmer"]]
        event_d["event_len"] = float(event_l[idx["event_len"]])
        if "start_idx" in idx:
            event_d["start_idx"] = int(event_l[idx["start_idx"]])
            event_d["end_idx"] = int(event_l[idx["end_idx"]])
        if "samples" in idx:
            event_d["sample_list"] = event_l[idx["samples"]].split(",")
        return event_d

    def _init_kmer_dict (self, event_d):
        """Start a new kmer dict from first event values"""
        kmer_d = OrderedDict ()
        kmer_d["ref_pos"] = event_d["ref_pos"]
        kmer_d["ref_kmer"] = event_d["ref_kmer"]
        kmer_d["num_events"] = 1
        kmer_d["dwell_time"] = event_d["event_len"]
        kmer_d["NNNNN_dwell_time"] = 0.0
        kmer_d["mismatch_dwell_time"] = 0.0
        if event_d["mod_kmer"] == "NNNNN":
            kmer_d["NNNNN_dwell_time"] += event_d["event_len"]
        elif event_d["mod_kmer"] != event_d["ref_kmer"]:
            kmer_d["mismatch_dwell_time"] += event_d["event_len"]
        if "start_idx" in event_d:
            kmer_d["start_idx"] = event_d["start_idx"]
            kmer_d["end_idx"] = event_d["end_idx"]
        if "sample_list" in event_d:
            kmer_d["sample_list"] = event_d["sample_list"]
        return kmer_d

    def _update_kmer_dict (self, kmer_d, event_d):
        """Update kmer dict from subsequent event values"""
        kmer_d["num_events"] += 1
        kmer_d["dwell_time"] += event_d["event_len"]
        if event_d["mod_kmer"] == "NNNNN":
            kmer_d["NNNNN_dwell_time"] += event_d["event_len"]
        elif event_d["mod_kmer"] != event_d["ref_kmer"]:
            kmer_d["mismatch_dwell_time"] += event_d["event_len"]
        if "start_idx" in event_d:
            kmer_d["start_idx"] = event_d["start_idx"]
        if "sample_list" in event_d:
            kmer_d["sample_list"].extend(event_d["sample_list"])
        return kmer_d

    def _kmer_dict_to_str (self, kmer_d):
        """"""
        # Write base fields
        s = "{}\t{}\t{}\t{}\t{}\t{}".format(
            kmer_d["ref_pos"],
            kmer_d["ref_kmer"],
            kmer_d["num_events"],
            kmer_d["dwell_time"],
            kmer_d["NNNNN_dwell_time"],
            kmer_d["mismatch_dwell_time"])
        # Facultative index fields
        if "start_idx" in kmer_d:
            s += "\t{}\t{}".format(
                kmer_d["start_idx"],
                kmer_d["end_idx"])
        # Facultative samples fields
        if "sample_list" in kmer_d:
            sample_array = np.array (kmer_d["sample_list"], dtype=np.float32)
            if "mean" in self.stat_fields:
                s += "\t{}".format(np.mean (sample_array))
            if "std" in self.stat_fields:
                s += "\t{}".format(np.std (sample_array))
            if "median" in self.stat_fields:
                s += "\t{}".format(np.median (sample_array))
            if "mad" in self.stat_fields:
                s += "\t{}".format(np.median(np.abs(sample_array-np.median(sample_array))))
            if "num_signals" in self.stat_fields:
                s += "\t{}".format(len(sample_array))
            if self.write_samples:
                s += "\t{}".format(",".join(kmer_d["sample_list"]))
        return s

    def _make_ouput_header (self, event_d):
        """"""
        # Write base fields
        s = "ref_pos\tref_kmer\tnum_events\tdwell_time\tNNNNN_dwell_time\tmismatch_dwell_time"
        # Write extra fields
        if "start_idx" in event_d:
            s += "\tstart_idx\tend_idx"
        if "sample_list" in event_d:
            if "mean" in self.stat_fields:
                s += "\tmean"
            if "std" in self.stat_fields:
                s += "\tstd"
            if "median" in self.stat_fields:
                s += "\tmedian"
            if "mad" in self.stat_fields:
                s += "\tmad"
            if "num_signals" in self.stat_fields:
                s += "\tnum_signals"
            if self.write_samples:
                s += "\tsamples"
        return s
