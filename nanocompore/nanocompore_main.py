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

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main ():
    # Simple triage function
    try:
        args = sys.argv
        if len(args) == 1:
            raise ValueError ("Error: Missing command\n")
        if args[1] == "sample_compare":
            sample_compare_main ()
        if args[1] == "model_compare":
            model_compare_main ()
        elif args[1] in ["-v", "--version"]:
            stderr_print ("{} v{}\n".format(package_name, package_version))
        elif args[1] in ["-h", "--help"]:
            raise ValueError ("nanocompore help\n")
        else:
            raise ValueError ("Error: Invalid command '{}'\n".format(args[1]))

    except ValueError as E:
        stderr_print (E)
        stderr_print ("Usage: nanocompore [command] [options]\n")
        stderr_print ("Valid command:\n\t-v/--version\n\tsample_compare\n\tmodel_compare\n")
        stderr_print ("For help on given command, type nanocompore [command] -h\n")
        sys.exit()

#~~~~~~~~~~~~~~SAMPLE COMPARE~~~~~~~~~~~~~~#
def sample_compare_main ():
    # Define parser object
    parser = argparse.ArgumentParser (description="Find differences in two nanopolish events files")
    parser.prog = "nanocompore sample_compare"

    # Define arguments
    parser.add_argument("subprogram")
    parser.add_argument("--file1", required=True,
        help="Path to the Nanopolish events file for sample 1")
    parser.add_argument("--file2", required=True,
        help="Path to the Nanopolish events file for sample 2")
    parser.add_argument("-o", "--outfolder", required=True,
        help="Path to a directory where to store the results. Must not exist")
    parser.add_argument("-n", type=int, default=6, required=False,
        help="Number of threads (two are used for reading files, all the other for processing in parallel).")
    parser.add_argument("--pthr", type=float, default=0.1, required=False,
        help="Adjusted p-value threshold for reporting sites.")
    parser.add_argument("--bedfile", default=None, required=False,
        help="Bedfile for annotating results")
    parser.add_argument("--min_coverage", type=int, default=10, required=False,
        help="Minimum number of reads covering a transcript (in each sample) for it to be considered. Set to 0 to consider all transcripts.")
    parser.add_argument("--logLevel", default="warning",
        help="Set the log level. Valid values: warning, info, debug.")
    # Parse Arguments
    a = parser.parse_args()

    # w = Whitelist ()
    # s = SampComp ()
    # s.write_results ()
    # s.do_other_stuff ()

#~~~~~~~~~~~~~~MODEL COMPARE~~~~~~~~~~~~~~#
def model_compare_main ():
    print ("Not implemented yet")
