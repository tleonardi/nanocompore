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
from nanocompore.PostProcess import PostProcess
from nanocompore.common import *

#~~~~~~~~~~~~~~MAIN PARSER ENTRY POINT~~~~~~~~~~~~~~#

def main():
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v' + package_version)
    subparsers = parser.add_subparsers(dest='subcommand',
        description=textwrap.dedent("""
            nanocompore implements the following subcommands\n
            \t* eventalign_collapse : Collapse the nanopolish eventalign output at kmer level and compute kmer-level statistics\n
            \t* sampcomp : Compare samples from two conditions and find significant signal differences\n
            \t* simreads : Simulate reads as a NanopolishComp-like file from a FASTA file and a built-in model"""))
    subparsers.required = True

    # Eventalign_collapse subparser
    parser_ec = subparsers.add_parser("eventalign_collapse", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        Collapse the nanopolish eventalign output at kmer level and compute kmer-level statistics
        * Minimal example:
            nanocompore eventalign_collapse -i nanopolish_eventalign.tsv -s T1"""))
    parser_ec.set_defaults(func=eventalign_collapse_main)

    parser_ec_in = parser_ec.add_argument_group("Input options")
    parser_ec_in.add_argument("--input1", "-1", required=True, nargs="+",
                              help="Filenames or paths of nanopolish eventalign TSV files for condition 1")
    parser_ec_in.add_argument("--input2", "-2", required=True, nargs="+",
                              help="Filenames or paths of nanopolish eventalign TSV files for condition 2")
    parser_ec_in.add_argument("--samples1", "-s1", nargs="*",
                              help="Sample names for input files from condition 1 (default: 'C1', 'C2', etc.)")
    parser_ec_in.add_argument("--samples2", "-s2", nargs="*",
                              help="Sample names for input files from condition 1 (default: 'T1', 'T2', etc.)")
    parser_ec_in.add_argument("--condition1", "-c1", default="control", help="Label for condition 1")
    parser_ec_in.add_argument("--condition2", "-c2", default="treatment", help="Label for condition 2")
    parser_ec_in.add_argument("--indir", "-i",
                              help="Path to input directory (will be preprended to 'input1' and 'input2')")
    # TODO: expose 'Eventalign_collapse' parameter 'master_db' here?

    parser_ec_out = parser_ec.add_argument_group("Output options")
    parser_ec_out.add_argument("--outdir", "-o", required=True, help="Directory for output files and subdirectories")
    parser_ec_out.add_argument("--overwrite", "-w", action="store_true",
                               help="Overwrite existing output files? (default: %(default)s)")
    parser_ec_out.add_argument("--subdirs", "-d", default=100, type=int,
                               help="Number of subdirectories to split output files into (default: %(default)s)")

    parser_ec_run = parser_ec.add_argument_group("Run options")
    parser_ec_run.add_argument("--n_lines", "-n", default=None , type=int ,
                               help = "Number of lines to parse (default: no limit)")

    parser_ec_misc = parser_ec.add_argument_group("Other options")
    parser_ec_misc.add_argument("--nthreads", "-t", default=3, type=int,
                                help="Total number of threads (minimum: 3, default: %(default)s)")

    # SampComp subparser
    parser_sc = subparsers.add_parser('sampcomp', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
            Compare 2 samples and find significant signal differences\n
            * Minimal example:\n
                nanocompore sampcomp -i eventalign_collapse.db -1 C1,C2 -2 T1,T2 -f ref.fa"""))
    parser_sc.set_defaults(func=sampcomp_main)

    parser_sc_in = parser_sc.add_argument_group('Input options')
    parser_sc_in.add_argument("--input", "-i", required=True,
        help="Path to the input/output directory containing 'eventalign_collapse' results (required)")
    parser_sc_in.add_argument("--fasta", "-f", required=True,
                              help="Fasta file used for mapping (required)")
    parser_sc_in.add_argument("--bed", default=None,
                              help="BED file with annotation of transcriptome used for mapping (optional)")

    parser_sc_out = parser_sc.add_argument_group("Output options")
    parser_sc_out.add_argument("--report", "-r", default="sampcomp.tsv",
                               help="Path or filename of report output file (default: %(default)s)")
    parser_sc_out.add_argument("--details", "-d", action="store_true",
                               help="Include additional details in report? (default: %(default)s)")

    parser_sc_filter = parser_sc.add_argument_group("Transcript filtering options")
    parser_sc_filter.add_argument("--max_invalid_kmers_freq", type=float, default=0.1,
        help="Maximum fequency of invalid kmers (default: %(default)s)")
    parser_sc_filter.add_argument("--min_coverage", type=int, default=30,
        help="Minimum coverage required in each condition to perform the comparison (default: %(default)s)")
    parser_sc_filter.add_argument("--downsample_high_coverage", type=int, default=5000,
        help="Downsample transcripts with high coverage to this number of reads (default: %(default)s)")
    parser_sc_filter.add_argument("--min_ref_length", type=int, default=100,
        help="Minimum length of a reference transcript for inclusion in the analysis (default: %(default)s)")

    parser_sc_test = parser_sc.add_argument_group('Statistical testing options')
    parser_sc_test.add_argument("--pvalue_threshold", "-p", type=float, default=0.01,
        help="Adjusted p-value threshold for reporting significant sites (default: %(default)s)")
    parser_sc_test.add_argument("--univariate_test", choices=["KS", "MW", "ST", "none"], default="KS",
        help="Univariate test for comparing kmer data between conditions. KS: Kolmogorov-Smirnov test, MW: Mann-Whitney test, ST: Student's t-test, none: no univariate test. (default: %(default)s)")
    parser_sc_test.add_argument("--no_gmm", action="store_true",
        help="Do not perform the GMM fit and subsequent test (see --gmm_test) (default: %(default)s)")
    parser_sc_test.add_argument("--gmm_test", choices=["logit", "anova", "none"], default="logit",
        help="Statistical test performed after GMM fitting (unless --no_gmm is used). (default: %(default)s)")
    parser_sc_test.add_argument("--store_gmm_components", choices=["all", "best", "none"], default="none",
        help="Store components of any fitted GMMs in the database? (default: %(default)s)")
    parser_sc_test.add_argument("--allow_warnings", action="store_true",
                                help="If True runtime warnings during the ANOVA tests (see --gmm_test) don't raise an error (default: %(default)s)")
    parser_sc_test.add_argument("--sequence_context", type=int, default=0, choices=range(0,5),
        help="Sequence context for combining p-values (default: %(default)s)")
    parser_sc_test.add_argument("--sequence_context_weights", default="uniform", choices=["uniform", "harmonic"],
        help="Type of position weighting to use for combining p-values (default: %(default)s)")

    parser_sc_misc = parser_sc.add_argument_group('Other options')
    parser_sc_misc.add_argument("--nthreads", "-t", type=int, default=3,
                                help="Number of threads (default: %(default)s)")

    # SimReads subparser
    parser_sr = subparsers.add_parser('simreads', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
            Simulate reads as a NanopolishComp like file from a fasta file and an inbuild model\n
            * Minimal example without model alteration
                nanocompore simreads -f ref.fa -o results -n 50\n
            * Minimal example with alteration of model intensity loc parameter for 50% of the reads
            nanocompore simreads -f ref.fa -o results -n 50 --intensity_mod 2 --mod_reads_freq 0.5 --mod_bases_freq 0.2"""))
    parser_sr.set_defaults(func=simreads_main)

    parser_sr_in = parser_sr.add_argument_group('Input options')
    parser_sr_in.add_argument("--fasta", "-f", required=True,
        help="FASTA file containing transcript sequences to use for artificial reads")

    parser_sr_out = parser_sr.add_argument_group("Output options")
    parser_sr_out.add_argument("--output", "-o", default="out",
                               help="Prefix for output files (default: %(default)s)")
    parser_sr_out.add_argument("--outdir", "-d", default="",
        help="Directory for output files. Will be preprended to --output if given. (default: %(default)s)")
    parser_sr_out.add_argument("--overwrite", "-w", action="store_true",
                               help="Overwrite existing output files? (default: %(default)s)")

    parser_sr_modify = parser_sr.add_argument_group('Signal modification options')
    parser_sr_modify.add_argument("--intensity_mod", type=float, default=0,
        help="Fraction of intensity distribution SD by which to modify the intensity distribution loc value (default: %(default)s)")
    parser_sr_modify.add_argument("--dwell_mod", type=float, default=0,
        help="Fraction of dwell time distribution SD by which to modify the intensity distribution loc value (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_reads_freq", type=float, default=0,
        help="Frequency of reads to modify (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_freq", type=float, default=0.25,
        help="Frequency of bases to modify in each read (if possible) (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_type", default="A", choices=["A","T","C","G"],
        help="Base for which to modify the signal (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_extend_context", type=int, default=2,
        help="number of adjacent base affected by the signal modification following an harmonic series (default: %(default)s)")
    parser_sr_modify.add_argument("--min_mod_dist", type=int, default=6,
        help="Minimal distance between 2 bases to modify (default: %(default)s)")
    parser_sr_misc = parser_sr.add_argument_group('Other options')
    parser_sr_misc.add_argument("--run_type", default="RNA", choices=["RNA", "DNA"],
        help="Define the run type model to import (default: %(default)s)")
    parser_sr_misc.add_argument("--nreads_per_ref", "-n", type=int, default=100,
        help="Number of reads to generate per references (default: %(default)s)")
    parser_sr_misc.add_argument("--pos_rand_seed", type=int, default=42 ,
        help="Define a seed for random position picking to get a deterministic behaviour (default: %(default)s)")
    parser_sr_misc.add_argument("--not_bound", action='store_true', default=False,
        help="Do not bind the values generated by the distributions to the observed min and max observed values from the model file (default: %(default)s)")

    # Add common options for all parsers
    for sp in [parser_ec, parser_sc, parser_sr]:
        sp_verbosity = sp.add_argument_group("Verbosity options")
        sp_verbosity.add_argument("--log_level", default="info", choices=["warning", "info", "debug", "trace"],
                                  help="Set the log level (default: %(default)s)")
    # only for SimReads:
    sp_verbosity.add_argument("--progress", action="store_true",
                              help="Display a progress bar during execution (default: %(default)s)")

    args = parser.parse_args()

    # Check if output folder already exists
    try:
        mkdir(fn=args.outdir, exist_ok=True)
    except (NanocomporeError, FileExistsError) as E:
        raise NanocomporeError(f"Could not create the output folder: {args.outdir}")

    # Set logger
    log_fn = os.path.join(args.outdir, vars(args)["subcommand"] + ".log")
    set_logger(args.log_level, log_fn=log_fn)

    # Call relevant subfunction
    args.func(args)

#~~~~~~~~~~~~~~SUBCOMMAND FUNCTIONS~~~~~~~~~~~~~~#

def eventalign_collapse_main(args):
    """"""
    logger.warning("Running Eventalign_collapse")

    # Init Eventalign_collapse
    if not args.samples1:
        args.samples1 = [f"C{n}" for n in range(1, len(args.input1) + 1)]
    elif len(args.samples1) != len(args.input1):
        raise NanocomporeError("Same number of arguments required for 'input1' and 'samples1'")
    if not args.samples2:
        args.samples2 = [f"T{n}" for n in range(1, len(args.input2) + 1)]
    elif len(args.samples2) != len(args.input2):
        raise NanocomporeError("Same number of arguments required for 'input2' and 'samples2'")
    if args.indir:
        if not os.path.isdir(indir):
            raise NanocomporeError(f"Input directory not found: {args.indir}")
        args.input1 = [os.path.join(indir, f) for f in args.input1]
        args.input2 = [os.path.join(indir, f) for f in args.input2]

    sample_dict = {args.condition1: [(f, s) for f, s in zip(args.input1, args.samples1)],
                   args.condition2: [(f, s) for f, s in zip(args.input2, args.samples2)]}

    e = Eventalign_collapse(sample_dict,
                            output_dir = args.outdir,
                            output_subdirs = args.subdirs,
                            overwrite = args.overwrite,
                            n_lines = args.n_lines,
                            nthreads = args.nthreads)

    # Run eventalign_collapse
    e()

def sampcomp_main(args):
    """"""
    logger.warning("Running SampComp")

    univar_test = args.univariate_test if args.univariate_test != "none" else None
    gmm_test = args.gmm_test if args.gmm_test != "none" else None
    thresholds = {}
    suffix = "_context" if args.sequence_context > 0 else ""
    if gmm_test:
        thresholds["gmm_pvalue" + suffix] = args.pvalue_threshold
    if univar_test:
        thresholds["intensity_pvalue" + suffix] = args.pvalue_threshold
        thresholds["dwell_pvalue" + suffix] = args.pvalue_threshold

    # Init SampComp
    s = SampComp(input_dir = args.input,
                 fasta_fn = args.fasta,
                 univariate_test = univar_test,
                 fit_gmm = not args.no_gmm,
                 gmm_test = gmm_test,
                 store_gmm_components = args.store_gmm_components,
                 allow_anova_warnings = args.allow_warnings,
                 sequence_context = args.sequence_context,
                 sequence_context_weighting = args.sequence_context_weights,
                 min_coverage = args.min_coverage,
                 min_transcript_length = args.min_ref_length,
                 downsample_high_coverage = args.downsample_high_coverage,
                 max_invalid_kmers_freq = args.max_invalid_kmers_freq,
                 significance_thresholds = thresholds,
                 nthreads = args.nthreads)

    # Run SampComp
    s()

    # Save all reports
    if not args.report:
        return

    adj_thresholds = {"adj_" + k: v for k, v in thresholds.items()}
    p = PostProcess(args.input, bed_path=args.bed)
    # TODO: update "save_all()" and call that instead
    p.save_report(args.report, adj_thresholds, args.details)


def simreads_main(args):
    """"""
    logger.warning("Running SimReads")

    # Run SimReads
    SimReads(fasta_fn = args.fasta,
             outpath = args.outdir,
             outprefix = args.output,
             overwrite = args.overwrite,
             run_type = args.run_type,
             nreads_per_ref = args.nreads_per_ref,
             intensity_mod = args.intensity_mod,
             dwell_mod = args.dwell_mod,
             mod_reads_freq = args.mod_reads_freq,
             mod_bases_freq = args.mod_bases_freq,
             mod_bases_type = args.mod_bases_type,
             mod_extend_context = args.mod_extend_context,
             min_mod_dist = args.min_mod_dist,
             pos_rand_seed = args.pos_rand_seed,
             not_bound = args.not_bound,
             progress = args.progress)

#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    # execute only if run as a script
    main()
