# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import namedtuple, Counter, OrderedDict
import logging

# Third party
import numpy as np
from tqdm import tqdm

# Local package
from nanocompore.helper_lib import counter_to_str, access_file, mytqdm
from nanocompore.NanocomporeError import NanocomporeError

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class Whitelist (object):

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__ (self,
        s1_index_fn,
        s2_index_fn,
        fasta_index_fn = None,
        min_coverage = 10,
        max_NNNNN_kmers_freq = 0.2,
        max_mismatching_kmers_freq = 0.2,
        max_missing_kmers_freq = 0.2,
        logLevel="info"):
        """
        s1_index_fn: Path to sample 1 eventalign_collapse index file
        s2_index_fn: Path to sample 2 eventalign_collapse index file
        fasta_index_fn: Path to a fasta index corresponding to the reference used for read alignemnt (see samtools faidx)
        min_coverage: minimal coverage required
        max_NNNNN_kmers_freq: maximum frequency of NNNNN kmers in reads (1 to deactivate)
        max_mismatching_kmers_freq: maximum frequency of mismatching kmers in reads (1 to deactivate)
        max_missing_kmers_freq: maximum frequency of missing kmers in reads (1 to deactivate)
        """

        # Set logging level
        logger.setLevel (logLevel_dict.get (logLevel, logging.WARNING))
        logger.info ("Initialise and checks options")
        self.__logLevel = logLevel

        # Check files
        for fn in (s1_index_fn, s2_index_fn, fasta_index_fn):
            access_file (fn)
        self.__s1_index_fn = s1_index_fn
        self.__s2_index_fn = s2_index_fn
        self.__fasta_index_fn = fasta_index_fn

        # Save other args
        self.__min_coverage = min_coverage
        self.__max_NNNNN_kmers_freq = max_NNNNN_kmers_freq
        self.__max_mismatching_kmers_freq = max_mismatching_kmers_freq
        self.__max_missing_kmers_freq = max_missing_kmers_freq

        # Read fasta index to get reference length
        logger.info ("Read fasta index files")
        self.ref_len_dict = self._read_fasta_index ()

        # Create reference index for both files
        logger.info ("Read eventalign index files")
        ref_reads = self._read_eventalign_index ()

        # First filtering pass at transcript level
        logger.info ("Filter out references with low coverage")
        ref_reads = self._select_ref (ref_reads)

        # Second filtering pass at base level
        logger.info ("Compute coverage per reference and select intervals with high enough coverage")
        self.ref_interval_reads = self._select_intervals (ref_reads)

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

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _read_fasta_index (self):
        """Read a fasta index file and return a refid:ref length dictionary"""

        ref_len_dict = OrderedDict ()
        with open (self.__fasta_index_fn) as fp:
            for line in fp:
                ls = line.rstrip().split()
                if len(ls) != 5:
                    raise NanocomporeError ("Invalid fasta index file: {}".format(fn))
                ref_len_dict[ls[0]] = int(ls[1])
        logger.info ("\tTotal references: {}".format(len(ref_len_dict)))
        return ref_len_dict

    def _read_eventalign_index (self):
        """Read the 2 index files and sort by sample and ref_id in a multi level dict"""

        ref_reads = OrderedDict ()

        for lab, fn in ("S1", self.__s1_index_fn), ("S2", self.__s2_index_fn):
            with open (fn) as fp:

                # get field names from header
                header = fp.readline().rstrip().split()
                for field in ("ref_id","ref_start","ref_end","read_id", "kmers", "NNNNN_kmers", "mismatching_kmers", "missing_kmers", "byte_offset", "byte_len"):
                    if not field in header:
                        raise NanocomporeError ("Invalid eventalign_collapse index file: {}".format(fn))

                line_tuple = namedtuple("line_tuple", header)
                c = Counter ()

                for line in fp:
                    ls = line.rstrip().split()
                    lt = line_tuple (ls[0], int(ls[1]), int(ls[2]), ls[3], int(ls[4]), int(ls[5]) , int(ls[6]) , int(ls[7]) , int(ls[8]), int(ls[9]))
                    # filter out reads with high number of problematic kmers
                    if lt.NNNNN_kmers/lt.kmers > self.__max_NNNNN_kmers_freq:
                        c ["high NNNNN_kmers reads"] += 1
                    elif lt.mismatching_kmers/lt.kmers > self.__max_mismatching_kmers_freq:
                        c ["high mismatching_kmers reads"] += 1
                    elif lt.missing_kmers/lt.kmers > self.__max_missing_kmers_freq:
                        c ["high missing_kmers reads"] += 1
                    # Save valid reads
                    else:
                        if not lt.ref_id in ref_reads:
                            ref_reads[lt.ref_id] = OrderedDict ()
                        if not lab in ref_reads [lt.ref_id]:
                            ref_reads[lt.ref_id][lab] = []
                        ref_reads[lt.ref_id][lab].append (lt)
                        c ["valid reads"] += 1

                logger.debug ("\tSample {} {}".format(lab, counter_to_str(c)))

        logger.info ("\tReferences found in index: {}".format(len(ref_reads)))
        return ref_reads

    def _select_ref (self, ref_reads):
        """Select ref_id with a minimal coverage in both sample"""

        selected_ref_reads = OrderedDict ()
        c = Counter()

        for ref_id, sample_reads in ref_reads.items ():
            if len(sample_reads) == 2 and len (sample_reads["S1"]) >= self.__min_coverage and len (sample_reads["S2"]) >= self.__min_coverage:
                selected_ref_reads [ref_id] = sample_reads
                if self.__logLevel == "debug":
                    c["ref_id"] += 1
                    c["positions"] += self.ref_len_dict[ref_id]
                    c["S1_reads"] += len (sample_reads["S1"])
                    c["S2_reads"] += len (sample_reads["S2"])

        logger.debug (counter_to_str(c))
        logger.info ("\tReferences remaining after reference coverage filtering: {}".format(len(selected_ref_reads)))
        return selected_ref_reads

    def _select_intervals (self, ref_reads):
        """Iterate over all transcripts find valid coverage intervals and extract overlapping reads"""

        ref_interval_reads = OrderedDict ()
        c = Counter ()

        pbar = tqdm (total=len(ref_reads), unit=" References", disable=self.__logLevel == "warning")
        for ref_id, sample_reads in ref_reads.items():
            pbar.update()

            # Compute reference coverage in a single numpy array
            cov_array = np.zeros ((2, self.ref_len_dict[ref_id]))
            for sample_index, (sample_id, read_list) in enumerate (sample_reads.items()):
                for read in read_list:
                    cov_array [sample_index][np.arange(read.ref_start, read.ref_end)] += 1

            # Get coordinates of intervals with minimum coverage
            valid_cov = False
            valid_interval_list = []
            for pos, (cov1, cov2) in enumerate (cov_array.T):
                # If coverage insuficient
                if cov1 < self.__min_coverage or cov2 < self.__min_coverage:
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
                for sample_id in ("S1", "S2"):
                    valid_reads = []
                    for read in sample_reads[sample_id]:
                        for interval_start, interval_end in valid_interval_list:
                            if read.ref_end >= interval_start and read.ref_start <= interval_end:
                                valid_reads.append (read)
                                if self.__logLevel == "debug":
                                    c["{}_reads".format(sample_id)] += 1
                                break
                    ref_interval_reads [ref_id] [sample_id] = valid_reads

        pbar.close ()
        logger.debug (counter_to_str(c))
        logger.info ("\tReferences remaining after position coverage filtering: {}".format(len(ref_interval_reads)))

        return ref_interval_reads
