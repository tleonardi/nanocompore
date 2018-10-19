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
        index_fn_dict,
        fasta_fn,
        min_coverage = 10,
        downsample_high_coverage = False,
        max_invalid_kmers_freq = 0.2,
        max_NNNNN_freq = 0.1,
        max_mismatching_freq = 0.1,
        max_missing_freq = 0.1,
        select_ref_id = [],
        exclude_ref_id = [],
        logLevel="info"):
        """
        index_fn_dict:
        fasta_fn: Path to a fasta file corresponding to the reference used for read alignemnt
        min_coverage: minimal coverage required in both samples
        downsample_high_coverage: For reference with higher coverage, downsample by randomly selecting reads.
        max_invalid_kmers_freq:
            If None, then the max_NNNNN_freq, max_mismatching_freq, max_missing_freq agrs will be used instead
        max_NNNNN_freq: maximum frequency of NNNNN kmers in reads (1 to deactivate)
        max_mismatching_freq: maximum frequency of mismatching kmers in reads (1 to deactivate)
        max_missing_freq: maximum frequency of missing kmers in reads (1 to deactivate)
        logLevel: Set the log level. Valid values: warning, info, debug
        """

        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise Whitelist and checks options")
        self.__logLevel = logLevel

        # Check index files
        for fn in index_fn_dict.values():
            if not access_file (fn):
                raise NanocomporeError("Cannot access nanopolish index file {}".format(fn))
            if not file_header_contains (fn, field_names=("ref_id","ref_start","ref_end","read_id","kmers","NNNNN_kmers","mismatching_kmers","missing_kmers","byte_offset","byte_len")):
                raise NanocomporeError("The index file {} does not contain the require header fields".format(fn))
        self.__index_fn_dict = index_fn_dict

        # Check fasta file
        if not access_file (fasta_fn):
            raise NanocomporeError("Cannot access fasta file {}".format(fasta_fn))
        # Read fasta index to get reference length
        try:
            logger.info ("Index Fasta file")
            self.__fasta = Fasta (fasta_fn)
        except:
            raise NanocomporeError("The fasta reference file cannot be read")

        # Create reference index for both files
        logger.info ("Read eventalign index files")
        ref_reads = self.__read_eventalign_index (
            index_fn_dict = index_fn_dict,
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
    def n_samples (self):
        return len(self.__index_fn_dict)

    @property
    def sample_id_list (self):
        return list (self.__index_fn_dict.keys())

    @property
    def ref_id_list (self):
        return list (self.ref_interval_reads.keys())

    @property
    def interval_dict (self):
        d = defaultdict (list)
        for red_id, ref_dict in self.ref_interval_reads.items():
            for interval in ref_dict["interval_list"]:
                d[red_id].append(interval)
        return d

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def __read_eventalign_index (self,
        index_fn_dict,
        max_invalid_kmers_freq,
        max_NNNNN_freq,
        max_mismatching_freq,
        max_missing_freq,
        select_ref_id,
        exclude_ref_id):
        """Read the 2 index files and sort by sample and ref_id in a multi level dict"""

        ref_reads = OrderedDict ()

        for lab, fn in index_fn_dict.items():
            with open (fn) as fp:

                # Get field names from header
                header = fp.readline().rstrip().split()
                line_tuple = namedtuple("line_tuple", header)

                c = Counter ()
                for line in fp:
                    lt = line_tuple (*numeric_cast_list(line.rstrip().split()))

                    # Filter out ref_id if a ref_id list was provided
                    if select_ref_id and not lt.ref_id in select_ref_id:
                        c ["Invalid ref_id"] += 1
                        continue
                    elif exclude_ref_id and lt.ref_id in exclude_ref_id:
                        c ["Invalid ref_id"] += 1
                        continue

                    # Filter out reads with high number of invalid kmers
                    if max_invalid_kmers_freq:
                        if (lt.NNNNN_kmers + lt.mismatching_kmers + lt.missing_kmers) / lt.kmers > max_invalid_kmers_freq:
                            c ["high invalid kmers reads"] += 1
                            continue
                    else:
                        if lt.NNNNN_kmers/lt.kmers > max_NNNNN_freq:
                            c ["high NNNNN_kmers reads"] += 1
                            continue
                        elif lt.mismatching_kmers/lt.kmers > max_mismatching_freq:
                            c ["high mismatching_kmers reads"] += 1
                            continue
                        elif lt.missing_kmers/lt.kmers > max_missing_freq:
                            c ["high missing_kmers reads"] += 1
                            continue

                    # Save valid reads
                    if not lt.ref_id in ref_reads:
                        ref_reads[lt.ref_id] = OrderedDict ()
                    if not lab in ref_reads [lt.ref_id]:
                        ref_reads[lt.ref_id][lab] = []
                    ref_reads[lt.ref_id][lab].append (lt)
                    c ["valid reads"] += 1

                logger.debug ("\tSample {} {}".format(lab, counter_to_str(c)))

        logger.info ("\tReferences found in index: {}".format(len(ref_reads)))
        return ref_reads

    def __select_ref (self, ref_reads, min_coverage):
        """Select ref_id with a minimal coverage in both sample"""

        valid_ref_reads = OrderedDict ()
        c = Counter()
        for ref_id, sample_reads in ref_reads.items ():
            try:
                # Verify that coverage is sufficient in all samples
                assert len(sample_reads) == self.n_samples
                for read_list in sample_reads.values():
                    assert len(read_list) >= min_coverage

                # Add to valid ref_id
                valid_ref_reads [ref_id] = sample_reads

                # Update Counter dict if debug mode
                if self.__logLevel == "debug":
                    c["valid_ref_id"] += 1
                    c["positions"] += len(self.__fasta[ref_id])
                    for sample_id, read_list in sample_reads.items():
                        c[sample_id+"_reads"] += len (read_list)

            except AssertionError:
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
        s_idx = {id:idx for idx, id in enumerate(self.sample_id_list)}
        ref_interval_reads = OrderedDict ()

        for ref_id, sample_reads in ref_reads.items():
            pbar.update()

            # Compute reference coverage in a single numpy array
            cov_array = np.zeros (shape=(len(self.__fasta[ref_id]), self.n_samples))
            for sample_id, read_list in sample_reads.items():
                sample_idx = s_idx [sample_id]
                for read in read_list:
                    cov_array [read.ref_start:read.ref_end+1, sample_idx] += 1

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

            # Intersect reads with valid coverage for both samples
            if valid_interval_list:
                if self.__logLevel == "debug":
                    c["ref_id"] += 1

                ref_interval_reads [ref_id] = OrderedDict ()
                ref_interval_reads [ref_id] ["interval_list"] = valid_interval_list

                # Select reads overlappig at least one valid interval
                for sample_id in self.sample_id_list:
                    valid_reads = []
                    for read in sample_reads [sample_id]:
                        for interval_start, interval_end in valid_interval_list:
                            if read.ref_end >= interval_start and read.ref_start <= interval_end:
                                valid_reads.append (read)
                                break

                    # Down sample if coverage too high
                    if downsample_high_coverage and len(valid_reads) > downsample_high_coverage:
                        valid_reads = random.sample (valid_reads, downsample_high_coverage)

                    if self.__logLevel == "debug":
                        c["{}_reads".format(sample_id)] += len(valid_reads)
                    ref_interval_reads [ref_id] [sample_id] = valid_reads

        pbar.close ()
        logger.debug (counter_to_str(c))
        logger.info ("\tReferences remaining after position coverage filtering: {}".format(len(ref_interval_reads)))

        return ref_interval_reads
