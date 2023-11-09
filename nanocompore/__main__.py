#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import argparse
from collections import *
import textwrap
import os

# Third party
from loguru import logger

# Local imports
from nanocompore import __version__ as package_version
from nanocompore import __name__ as package_name
from nanocompore import __description__ as package_description
from nanocompore.SampComp import SampComp
from nanocompore.SimReads import SimReads
from nanocompore.Eventalign_collapse import Eventalign_collapse
from nanocompore.common import *

#~~~~~~~~~~~~~~MAIN PARSER ENTRY POINT~~~~~~~~~~~~~~#

def main(args=None):
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v'+package_version)
    subparsers = parser.add_subparsers(dest='subcommand',
        description=textwrap.dedent("""
            nanocompore implements the following subcommands\n
            \t* sampcomp : Compare 2 samples and find significant signal differences\n"""))
    subparsers.required = True

    # Sampcomp subparser
    parser_sc = subparsers.add_parser('sampcomp', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
            Compare 2 samples and find significant signal differences\n
            * Minimal example \n
                nanocompore sampcomp -s samples.tsv -f ref -o results"""))
    parser_sc.set_defaults(func=sampcomp_main)
    parser_sc.add_argument("--sample_tsv", "-s", default=None, type=str, metavar="sample_tsv",
        help="TSV file containing the sample file labels. One sample per line with tabs separating " 
             "Sample_label, Condition_label, path/to/pod5, and path/to/bam. See full formatting description in documentation.")

    parser_sc_io = parser_sc.add_argument_group('Input options')
    parser_sc_io.add_argument("--fasta", "-f", type=str, required=True,
        help="Reference fasta file that was used for mapping the data (required)")
    parser_sc_io.add_argument("--bed", type=str, default=None,
        help="BED file with annotation of transcriptome used for mapping (optional)")

    parser_sc_filtering = parser_sc.add_argument_group('Transcript filtering options')
    parser_sc_filtering.add_argument("--min_coverage", type=int, default=30,
        help="Minimum coverage required in each condition to do the comparison (default: %(default)s)")
    parser_sc_filtering.add_argument("--downsample_high_coverage", type=int, default=5000,
        help="Transcripts with high coverage will be downsampled (default: %(default)s)")
    parser_sc_filtering.add_argument("--min_ref_length", type=int, default=100,
        help="Minimum length of a reference transcript to include it in the analysis (default: %(default)s)")

    parser_sc_testing = parser_sc.add_argument_group('Statistical testing options')
    parser_sc_testing.add_argument("--comparison_methods", type=str, default="GMM,KS",
        help="Comma separated list of comparison methods. Valid methods are: GMM,KS,TT,MW. (default: %(default)s)")
    parser_sc_testing.add_argument("--sequence_context", type=int, default=0, choices=range(0,5),
        help="Sequence context for combining p-values (default: %(default)s)")
    parser_sc_testing.add_argument("--sequence_context_weights", type=str, default="uniform", choices=["uniform", "harmonic"],
        help="Type of weights to use for combining p-values")
    parser_sc_testing.add_argument("--pvalue_thr", type=float, default=0.05,
        help="Adjusted p-value threshold for reporting significant sites (default: %(default)s)")
    parser_sc_testing.add_argument("--logit", action='store_true',
        help="Use logistic regression testing downstream of GMM method. This is a legacy option and is now the deault.")
    parser_sc_testing.add_argument("--anova", action='store_true',
        help="Use Anova test downstream of GMM method (default: %(default)s)")
    parser_sc_testing.add_argument("--allow_warnings", action='store_true', default=False,
        help="If True runtime warnings during the ANOVA tests don't raise an error (default: %(default)s)")

    parser_sc_misc = parser_sc.add_argument_group('Other options')
    parser_sc_misc.add_argument("--nthreads", "-t", type=int, default=3,
        help="Number of threads (default: %(default)s)")

    # Add common options for all parsers
    for sp in [parser_sc]:
        sp_output = sp.add_argument_group("Output options")
        sp_output.add_argument("--outpath", "-o", type=str, default="./",
            help="Path to the output folder (default: %(default)s)")
        sp_output.add_argument("--outprefix", "-p", type=str, default="out",
            help="text outprefix for all the files generated (default: %(default)s)")
        sp_output.add_argument("--overwrite", "-w", action='store_true', default=False,
            help="Use --outpath even if it exists already (default: %(default)s)")
        sp_verbosity = sp.add_argument_group("Verbosity options")
        sp_verbosity.add_argument("--log_level", type=str, default="info", choices=["warning", "info", "debug"],
            help="Set the log level (default: %(default)s)")
        sp_verbosity.add_argument("--progress", default=False, action='store_true',
            help="Display a progress bar during execution (default: %(default)s)")

    # Parse agrs and
    args = parser.parse_args()

    # Check if output folder already exists
    try:
        mkdir(fn=args.outpath, exist_ok=args.overwrite)
    except (NanocomporeError, FileExistsError) as E:
        raise NanocomporeError("Could not create the output folder. Try using `--overwrite` option or use another directory")

    # Set logger
    log_fn = os.path.join(args.outpath, args.outprefix+"_{}.log".format(vars(args)["subcommand"]))
    set_logger(args.log_level, log_fn=log_fn)

    # Call relevant subfunction
    args.func(args)

#~~~~~~~~~~~~~~SUBCOMMAND FUNCTIONS~~~~~~~~~~~~~~#

def sampcomp_main(args):
    """"""
    logger.warning("Running SampComp")

    # Init SampComp
    s = SampComp(
        in_tsv = args.sample_tsv,
        outpath = args.outpath,
        outprefix = args.outprefix,
        overwrite = args.overwrite,
        fasta_fn = args.fasta,
        bed_fn = args.bed,
        nthreads = args.nthreads,
        min_coverage = args.min_coverage,
        min_ref_length = args.min_ref_length,
        downsample_high_coverage = args.downsample_high_coverage,
        comparison_methods = args.comparison_methods,
        logit = True,
        anova = args.anova,
        allow_warnings = args.allow_warnings,
        sequence_context = args.sequence_context,
        sequence_context_weights = args.sequence_context_weights,
        progress = args.progress)

    # Run SampComp
    s()
#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    # execute only if run as a script
    main()
