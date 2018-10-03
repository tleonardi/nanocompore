#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
from nanocompore import __version__ as package_version
from nanocompore import __name__ as package_name
from nanocompore.SampComp import SampComp
from nanocompore.Whitelist import Whitelist
from nanocompore.common import NanocomporeError


#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main ():
    # Simple triage function
    try:
        args = sys.argv
        if len(args) == 1:
            raise NanocomporeError ("Error: Missing command\n")
        elif args[1] == "sample_compare":
            sample_compare_main ()
        elif args[1] == "model_compare":
            model_compare_main ()
        elif args[1] in ["-v", "--version"]:
            print ("{} v{}\n".format(package_name, package_version))
        elif args[1] in ["-h", "--help"]:
            raise NanocomporeError ("nanocompore help\n")
        else:
            raise NanocomporeError ("Error: Invalid command '{}'\n".format(args[1]))

    except NanocomporeError as E:
        print (E)
        print ("Usage: nanocompore [command] [options]\n")
        print ("Valid command:\n\t-v/--version\n\tsample_compare\n\tmodel_compare\n")
        print ("For help on given command, type nanocompore [command] -h\n")
        sys.exit()

#~~~~~~~~~~~~~~SAMPLE COMPARE~~~~~~~~~~~~~~#
def sample_compare_main ():
    # Define parser object
    parser = argparse.ArgumentParser (description="Find differences in two nanopolish eventalign collapsed files")
    parser.prog = "nanocompore sample_compare"

    # Define arguments
    parser.add_argument("subprogram")

    parser.add_argument("-1", "--s1_fn", required=True,
        help= "Path to sample 1 eventalign_collapse data file")
    parser.add_argument("-2", "--s2_fn", required=True,
        help= "Path to sample 2 eventalign_collapse data file")
    parser.add_argument("-f", "--fasta_index_fn", required=True,
        help= "Path to a fasta index corresponding to the reference used for read alignemnt (see samtools faidx)")
    parser.add_argument("-o", "--output_db_fn", required=True,
        help= "Path where to write the result database")
    parser.add_argument("--min_coverage", default=10 , type=int,
        help= "minimal coverage required in both samples")
    parser.add_argument("--downsample_high_coverage", default=0 , type=int,
        help= "For reference with higher coverage, downsample by randomly selecting reads.")
    parser.add_argument("--max_NNNNN_freq", default=0.2 , type=float,
        help= "maximum frequency of NNNNN kmers in reads (1 to deactivate)")
    parser.add_argument("--max_mismatching_freq", default=0.2 , type=float,
        help= "maximum frequency of mismatching kmers in reads (1 to deactivate)")
    parser.add_argument("--max_missing_freq", default=0.2 , type=float,
        help= "maximum frequency of missing kmers in reads (1 to deactivate)")
    parser.add_argument("--padj_threshold", default=0.1 , type=float,
        help= "Adjusted p-value threshold for reporting sites")
    parser.add_argument("--comparison_method", default="kmean",
        help= "Statistical method to compare the 2 samples signal")
    parser.add_argument("--sequence_context", default=0 , type=int,
        help= "Extend statistical analysis to contigous adjacent base is available")
    parser.add_argument("-t", "--nthreads", default=4 , type=int,
        help= "Number of threads, 2 are used for reading and writing, all the others for processing in parallel")
    parser.add_argument("--logLevel", default="info",
        help= "Set the log level. Valid values: warning, info, debug")

    a = parser.parse_args()

    s = SampComp(
        s1_fn = a.s1_fn,
        s2_fn = a.s2_fn,
        output_db_fn = a.output_db_fn,
        fasta_index_fn = a.fasta_index_fn,
        whitelist = a.whitelist,
        padj_threshold = a.padj_threshold,
        comparison_method = a.comparison_method,
        sequence_context = a.sequence_context,
        min_coverage = a.min_coverage,
        downsample_high_coverage = a.downsample_high_coverage,
        max_NNNNN_freq = a.max_NNNNN_freq,
        max_mismatching_freq = a.max_mismatching_freq,
        max_missing_freq = a.max_missing_freq,
        nthreads = a.nthreads,
        logLevel = a.logLevel)

    # s.write_results ()
    # s.do_other_stuff ()

#~~~~~~~~~~~~~~MODEL COMPARE~~~~~~~~~~~~~~#
def model_compare_main ():
    print ("Not implemented yet")
