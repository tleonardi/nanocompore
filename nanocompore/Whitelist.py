# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import namedtuple, Counter, OrderedDict, defaultdict
import logging
import random

# Third party
import numpy as np
from tqdm import tqdm
from pyfaidx import Fasta

# Local package
from nanocompore.common import counter_to_str, access_file, NanocomporeError, file_header_contains, numeric_cast_list

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist (object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__ (self,
        eventalign_fn_dict,
        fasta_fn,
        min_coverage = 10,
        downsample_high_coverage = False,
        max_invalid_kmers_freq = 0.1,
        max_NNNNN_freq = 0.1,
        max_mismatching_freq = 0.1,
        max_missing_freq = 0.1,
        select_ref_id = [],
        exclude_ref_id = [],
        logLevel="info"):
        """
        eventalign_fn_dict: Multilevel dictionnary indicating the condition_label, sample_label and file name of the eventalign_collapse output
            example d = {"S1": {"R1":"path1.tsv", "R2":"path2.tsv"}, "S2": {"R1":"path3.tsv", "R2":"path4.tsv"}}
        fasta_fn: Path to a fasta file corresponding to the reference used for read alignemnt
        min_coverage: minimal coverage required in both samples
        downsample_high_coverage: For reference with higher coverage, downsample by randomly selecting reads.
        max_invalid_kmers_freq: maximum frequency of NNNNN, mismatching and missing kmers in reads
            If None, then the max_NNNNN_freq, max_mismatching_freq, max_missing_freq agrs will be used instead
        max_NNNNN_freq: maximum frequency of NNNNN kmers in reads (1 to deactivate)
        max_mismatching_freq: maximum frequency of mismatching kmers in reads (1 to deactivate)
        max_missing_freq: maximum frequency of missing kmers in reads (1 to deactivate)
        select_ref_id: if given, only reference ids in the list will be selected for the analysis
        exclude_ref_id: if given, refid in the list will be excluded from the analysis
        logLevel: Set the log level. Valid values: warning, info, debug
        """

        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise Whitelist and checks options")
        self.__logLevel = logLevel

        # Check index files
        for sample_dict in eventalign_fn_dict.values():
            for fn in sample_dict.values():
                idx_fn = fn+".idx"
                if not access_file (idx_fn):
                    raise NanocomporeError("Cannot access eventalign_collapse index file {}".format(idx_fn))
                if not file_header_contains (idx_fn, field_names=("ref_id","ref_start","ref_end","read_id","kmers","NNNNN_kmers","mismatching_kmers","missing_kmers","byte_offset","byte_len")):
                    raise NanocomporeError("The index file {} does not contain the require header fields".format(idx_fn))
        self.__eventalign_fn_dict = eventalign_fn_dict

        # Get number of samples
        n = 0
        for sample_dict in self.__eventalign_fn_dict.values():
            for sample_lab in sample_dict.keys():
                n+=1
        self.__n_samples = n

        # Test is Fasta can be opened
        try:
            with Fasta (fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        # Create reference index for both files
        logger.info ("Read eventalign index files")
        ref_reads = self.__read_eventalign_index (
            eventalign_fn_dict = eventalign_fn_dict,
            max_invalid_kmers_freq = max_invalid_kmers_freq,
            max_NNNNN_freq = max_NNNNN_freq,
            max_mismatching_freq = max_mismatching_freq,
            max_missing_freq = max_missing_freq,
            select_ref_id = select_ref_id,
            exclude_ref_id = exclude_ref_id)

        # First filtering pass at transcript level
        logger.info ("Filter out references with low coverage")
        ref_reads = self.__select_ref (
            ref_reads = ref_reads,
            min_coverage=min_coverage)

        # Second filtering pass at base level
        logger.info ("Compute coverage per reference and select intervals with high enough coverage")
        self.ref_interval_reads = self.__select_intervals (
            ref_reads = ref_reads,
            min_coverage = min_coverage,
            downsample_high_coverage=downsample_high_coverage)

        self.ref_reads = ref_reads

    def __repr__ (self):
        m = "Whitelist: Number of references: {}".format(len(self))
        return m

    def __len__ (self):
        return len(self.ref_interval_reads)

    def __iter__ (self):
        for i, j in self.ref_interval_reads.items():
            yield (i,j)

    def __getitem__(self, items):
        return self.ref_interval_reads.get(items, None)

    #~~~~~~~~~~~~~~PUBLIC METHODS AND PROPERTIES~~~~~~~~~~~~~~#
    def to_bed (self, bed_fn):
        with open (bed_fn, "w") as fp:
            for ref_id, ref_dict in self.ref_interval_reads.items ():
                for start, end in ref_dict["interval_list"]:
                    fp.write ("{}\t{}\t{}\n".format (ref_id, start, end))

    @property
    def ref_id_list (self):
        return list (self.ref_interval_reads.keys())

    @property
    def interval_dict (self):
        d = defaultdict (list)
        for ref_id, ref_dict in self.ref_interval_reads.items():
            for interval in ref_dict["interval_list"]:
                d[ref_id].append(interval)
        return d

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def __read_eventalign_index (self,
        eventalign_fn_dict,
        max_invalid_kmers_freq,
        max_NNNNN_freq,
        max_mismatching_freq,
        max_missing_freq,
        select_ref_id,
        exclude_ref_id):
        """Read the 2 index files and sort by sample and ref_id in a multi level dict"""

        ref_reads = OrderedDict ()

        for cond_lab, sample_dict in eventalign_fn_dict.items():
            for sample_lab, fn in sample_dict.items ():
                idx_fn = fn+".idx"
                with open (idx_fn) as fp:

                    # Get column names from header
                    col_names = fp.readline().rstrip().split()
                    c = Counter ()
                    for line in fp:
                        try:
                            # Transform line to dict and cast str numbers to actual numbers
                            read = dict (zip(col_names, numeric_cast_list(line.rstrip().split())))

                            # Filter out ref_id if a select_ref_id list or exclude_ref_id list was provided
                            if select_ref_id and not read["ref_id"] in select_ref_id:
                                raise NanocomporeError ("Ref_id not in select list")
                            elif exclude_ref_id and read["ref_id"] in exclude_ref_id:
                                raise NanocomporeError ("Ref_id in exclude list")

                            # Filter out reads with high number of invalid kmers
                            if max_invalid_kmers_freq:
                                if(read["NNNNN_kmers"]+read["mismatching_kmers"]+read["missing_kmers"])/read["kmers"] > max_invalid_kmers_freq:
                                    raise NanocomporeError ("High invalid kmers reads")
                            else:
                                if max_NNNNN_freq and read["NNNNN_kmers"]/read["kmers"] > max_NNNNN_freq:
                                    raise NanocomporeError ("High NNNNN kmers reads")
                                elif max_mismatching_freq and read["mismatching_kmers"]/read["kmers"] > max_mismatching_freq:
                                    raise NanocomporeError ("High mismatching_kmers reads")
                                elif max_missing_freq and read["missing_kmers"]/read["kmers"] > max_missing_freq:
                                    raise NanocomporeError ("High missing_kmers reads")

                            # Create dict arborescence and save valid reads
                            if not read["ref_id"] in ref_reads:
                                ref_reads[read["ref_id"]] = OrderedDict()
                            if not cond_lab in ref_reads[read["ref_id"]]:
                                ref_reads[read["ref_id"]][cond_lab] = OrderedDict()
                            if not sample_lab in ref_reads[read["ref_id"]][cond_lab]:
                                ref_reads[read["ref_id"]][cond_lab][sample_lab] = []

                            # Fill in list of reads
                            ref_reads[read["ref_id"]][cond_lab][sample_lab].append (read)
                            c ["valid reads"] += 1

                        except NanocomporeError as E:
                            c [str(E)] += 1

                logger.debug ("\tCondition:{} Sample:{} {}".format(cond_lab, sample_lab, counter_to_str(c)))

        logger.info ("\tReferences found in index: {}".format(len(ref_reads)))
        return ref_reads

    def __select_ref (self,
        ref_reads,
        min_coverage):
        """Select ref_id with a minimal coverage in both sample"""

        valid_ref_reads = OrderedDict ()
        c = Counter()
        with Fasta (self._fasta_fn) as fasta:
            for ref_id, cond_dict in ref_reads.items():
                n=0
                for cond_lab, sample_dict in cond_dict.items():
                    for sample_lab, read_list in sample_dict.items():
                        if len(read_list) >= min_coverage:
                            n+=1

                # Add to valid ref_id
                if n == self.__n_samples:
                    valid_ref_reads [ref_id] = cond_dict
                    # Update Counter dict if debug mode
                    if self.__logLevel == "debug":
                        c["valid_ref_id"] += 1
                        c["positions"] += len(fasta[ref_id])
                        for cond_lab, sample_dict in cond_dict.items():
                            for sample_lab, read_list in sample_dict.items():
                                lab = "{} {} Reads".format(cond_lab, sample_lab)
                                c[lab] += len(read_list)
                else:
                    c["invalid_ref_id"] += 1

        logger.debug (counter_to_str(c))
        logger.info ("\tReferences remaining after reference coverage filtering: {}".format(len(valid_ref_reads)))
        return valid_ref_reads

    def __select_intervals (self,
         ref_reads,
         min_coverage,
         downsample_high_coverage):
        """Iterate over all transcripts find valid coverage intervals and extract overlapping reads"""

        c = Counter ()
        pbar = tqdm (total=len(ref_reads), unit=" References", disable=self.__logLevel == "warning")

        # Create a dict to associate each sample to an index in the coverage array
        idx = OrderedDict()
        n = 0
        for cond_lab, sample_dict in self.__eventalign_fn_dict.items():
            idx[cond_lab] = OrderedDict()
            for sample_lab in sample_dict.keys():
                idx[cond_lab][sample_lab] = n
                n+=1

        # Iterate over the dict cobtaining valid reads per conditions and samples
        ref_interval_reads = OrderedDict ()
        with Fasta (self._fasta_fn) as fasta:
            for ref_id, cond_dict in ref_reads.items():
                pbar.update()

                # Compute reference coverage in a single numpy array
                cov_array = np.zeros (shape=(len(fasta[ref_id]), self.__n_samples))
                for cond_lab, sample_dict in cond_dict.items():
                    for sample_lab, read_list in sample_dict.items():
                        sample_idx = idx[cond_lab][sample_lab]
                        for read in read_list:
                            cov_array [read["ref_start"]:read["ref_end"]+1, sample_idx] += 1

                # Reduce to min coverage per position
                cov_array = cov_array.min (axis=1)

                # Get coordinates of intervals with minimum coverage
                valid_cov = False
                valid_interval_list = []
                for pos, cov in enumerate (cov_array):
                    # If coverage insuficient
                    if cov < min_coverage:
                        if valid_cov:
                            valid_interval_list.append ((ref_start, ref_end))
                            if self.__logLevel == "debug":
                                c["intervals"] += 1
                                c["positions"] += abs (ref_start-ref_end)
                        valid_cov = False
                    # If the coverage is high enough for both samples
                    else:
                        if valid_cov:
                            ref_end = pos
                        else:
                            ref_start = ref_end = pos
                            valid_cov = True

                # Last valid interval exception
                if valid_cov:
                    valid_interval_list.append ((ref_start, ref_end))
                    if self.__logLevel == "debug":
                        c["intervals"] += 1
                        c["positions"] += abs (ref_start-ref_end)

                # Intersect reads with valid coverage for all samples
                if valid_interval_list:
                    if self.__logLevel == "debug":
                        c["ref_id"] += 1
                    ref_interval_reads [ref_id] = OrderedDict ()
                    ref_interval_reads [ref_id] ["interval_list"] = valid_interval_list

                    # Select reads overlapping at least one valid interval
                    for cond_lab, sample_dict in cond_dict.items():
                        for sample_lab, read_list in sample_dict.items():
                            valid_reads = []
                            for read in read_list:
                                for interval_start, interval_end in valid_interval_list:
                                    if read["ref_end"] >= interval_start and read["ref_start"] <= interval_end:
                                        valid_reads.append (read)
                                        break
                            # Downsample if coverage too high
                            if downsample_high_coverage and len(valid_reads) > downsample_high_coverage:
                                valid_reads = random.sample (valid_reads, downsample_high_coverage)

                            if self.__logLevel == "debug":
                                lab = "{} {} Reads".format(cond_lab, sample_lab)
                                c[lab] += len(valid_reads)

                            if not cond_lab in ref_interval_reads[ref_id]:
                                ref_interval_reads[ref_id][cond_lab] = OrderedDict ()
                            ref_interval_reads[ref_id][cond_lab][sample_lab] = valid_reads

        pbar.close ()
        logger.debug (counter_to_str(c))
        logger.info ("\tReferences remaining after position coverage filtering: {}".format(len(ref_interval_reads)))

        return ref_interval_reads
