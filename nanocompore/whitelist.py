# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import namedtuple, Counter, OrderedDict

# Third party
import numpy as np

#~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~#

def whitelist (
    s1_fn,
    s2_fn,
    fasta_index_fn = None,
    min_cov = 10,
    max_NNNNN_kmers_freq = 0.2,
    max_mismatching_kmers_freq = 0.2,
    max_missing_kmers_freq = 0.2,
    verbose=True):
    """
    * s1_fn
        Path to sample 1 eventalign_collapse index file
    * s2_fn
        Path to sample 2 eventalign_collapse index file
    * fasta_index_fn
        Path to a fasta index corresponding to the reference used for read alignemnt
        See Samtools faidx (http://www.htslib.org/doc/faidx.html)
    * min_cov
        minimal coverage required
    * max_NNNNN_kmers_freq
        maximum frequency of NNNNN kmers in reads
    * max_mismatching_kmers_freq
        maximum frequency of mismatching kmers in reads
    * max_missing_kmers_freq
        maximum frequency of missing kmers in reads
    """

    # Read fasta index
    if verbose: print ("Read fasta index")
    ref_len_dict = _read_fasta_index (fn=fasta_index_fn)
    if verbose: print (f"\tTotal references: {len(ref_len_dict)}")

    # Create reference index for both files
    if verbose: print ("Read eventalign index files")
    ref_reads = _read_eventalign_index (s1_fn, s2_fn, max_NNNNN_kmers_freq, max_mismatching_kmers_freq, max_missing_kmers_freq)
    if verbose: print (f"\tTotal references found {len(ref_reads)}")

    # Intersect both samples
    if verbose: print ("Filter out references with low coverage")
    ref_reads = _select_ref (ref_reads=ref_reads, min_cov=min_cov)
    if verbose: print (f"\tTranscripts remaining after reference coverage filtering: {len(ref_reads)}")

    if verbose: print ("Compute coverage per reference and select intervals with high enough coverage")
    ref_interval_reads = OrderedDict ()
    for ref_id, sample_reads in ref_reads.items ():
        # Compute reference coverage
        cov_array = _compute_ref_cov (sample_reads=sample_reads, ref_len=ref_len_dict[ref_id])
        # Get coordinates of intervals with minimum coverage
        valid_interval_list = _get_valid_intervals (cov_array, min_cov)
        # Intesect reads with valid coverage for both samples
        ref_interval_reads [ref_id] = _intersect_reads_interval (valid_interval_list, sample_reads)

    return ref_interval_reads

#~~~~~~~~~~~~~~HELPER PRIVATE FUNCTIONS~~~~~~~~~~~~~~#
def _read_fasta_index (fn):
    """Read a fasta index file and return a refid:ref length dictionary"""
    ref_len = OrderedDict ()
    with open (fn) as fp:
        for line in fp:
            ls = line.rstrip().split()
            ref_len[ls[0]] = int(ls[1])
    return ref_len

def _read_eventalign_index (s1_fn, s2_fn, max_NNNNN_kmers_freq, max_mismatching_kmers_freq, max_missing_kmers_freq):
    """Read the 2 index files and sort by sample and ref_id in a multi level dict"""
    ref_reads = OrderedDict ()

    for lab, fn in ("S1", s1_fn), ("S2", s2_fn):
        with open (fn) as fp:
            # get field names from header
            header = fp.readline().rstrip().split()
            line_tuple = namedtuple("line_tuple", header)
            c = Counter ()
            for line in fp:
                ls = line.rstrip().split()
                lt = line_tuple (ls[0], int(ls[1]), int(ls[2]), ls[3], int(ls[4]), int(ls[5]) , int(ls[6]) , int(ls[7]) , int(ls[8]))
                # filter out reads with high number of problematic kmers
                if max_NNNNN_kmers_freq and lt.NNNNN_kmers/lt.kmers > max_NNNNN_kmers_freq:
                    c ["high NNNNN_kmers reads"] += 1
                elif max_mismatching_kmers_freq and lt.mismatching_kmers/lt.kmers > max_mismatching_kmers_freq:
                    c ["high mismatching_kmers reads"] += 1
                elif max_missing_kmers_freq and lt.missing_kmers/lt.kmers > max_missing_kmers_freq:
                    c ["high missing_kmers reads"] += 1
                # Save valid reads
                else:
                    if not lt.ref_id in ref_reads:
                        ref_reads[lt.ref_id] = OrderedDict ()
                    if not lab in ref_reads [lt.ref_id]:
                        ref_reads[lt.ref_id][lab] = []
                    ref_reads[lt.ref_id][lab].append (lt)
                    c ["valid reads"] += 1
        print (c)
    return ref_reads

def _select_ref (ref_reads, min_cov):
    """Select ref_id with a minimal coverage in both sample"""
    invalid_ref = []
    for ref_id, sample_reads in ref_reads.items ():
        if len(sample_reads) < 2 or len (sample_reads["S1"]) < min_cov or len (sample_reads["S2"]) < min_cov:
            invalid_ref.append (ref_id)
    for ref_id in invalid_ref:
        del ref_reads [ref_id]
    return ref_reads

def _compute_ref_cov (sample_reads, ref_len):
    """Compute the base level coverage for a reference for both samples"""
    cov_array = np.zeros ((2, ref_len))
    for read in sample_reads["S1"]:
        cov_array [0][np.arange(read.ref_start, read.ref_end)] += 1
    for read in sample_reads["S2"]:
        cov_array [1][np.arange(read.ref_start, read.ref_end)] += 1
    return cov_array

def _get_valid_intervals (cov_array, min_cov):
    """Find intervals with a coverage higher than min_cov and return a list"""
    valid_cov = False
    valid_interval_list = []
    for pos, (cov1, cov2) in enumerate (cov_array.T):
        # If coverage insuficient
        if cov1 < min_cov or cov2 < min_cov:
            if valid_cov:
                valid_interval_list.append ((ref_start, ref_end))
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

    return valid_interval_list

def _intersect_reads_interval (valid_interval_list, sample_reads):
    """Get reads intesecting a particular interval and return a multi level dict"""
    ref_interval_reads = OrderedDict ()
    for interval_start, interval_end in valid_interval_list:
        ref_interval_reads [(interval_start, interval_end)] = {"S1":[], "S2":[]}

    for sample_id, read_list in sample_reads.items():
        for read in read_list:
            for interval_start, interval_end in valid_interval_list:
                if read.ref_end >= interval_start and read.ref_start <= interval_end:
                    ref_interval_reads[(interval_start, interval_end)][sample_id].append (read)
    return ref_interval_reads
