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
            \t* eventalign_collapse : Collapse the nanopolish eventalign output at kmers level and compute kmer level statistics\n
            \t* sampcomp : Compare 2 samples and find significant signal differences\n
            \t* simreads : Simulate reads as a NanopolishComp like file from a fasta file and an inbuild model"""))
    subparsers.required = True

    # Sampcomp subparser
    parser_sc = subparsers.add_parser('sampcomp', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
            Compare 2 samples and find significant signal differences\n
            * Minimal example with file_list arguments\n
                nanocompore sampcomp -1 f1.tsv,f2.tsv -2 f3.tsv,f4.tsv -f ref.fa -o results
            * Minimal example with sample YAML file\n
                nanocompore sampcomp -y samples.yaml -f ref -o results"""))
    parser_sc.set_defaults(func=sampcomp_main)
    parser_sc_sample_yaml = parser_sc.add_argument_group('YAML sample files', description="Option allowing to describe sample files in a YAML file")
    parser_sc_sample_yaml.add_argument("--sample_yaml", "-y", default=None, type=str, metavar="sample_yaml",
        help="YAML file containing the sample file labels. See formatting in documentation. (required if --file_list1 and --file_list2 not given)")
    parser_sc_sample_args = parser_sc.add_argument_group('Arguments sample files', description="Option allowing to describe sample files directly as command line arguments")
    parser_sc_sample_args.add_argument("--file_list1", "-1", default=None, type=str, metavar="/path/to/Condition1_rep1,/path/to/Condition1_rep2",
        help="Comma separated list of NanopolishComp files for label 1. (required if --sample_yaml not given)")
    parser_sc_sample_args.add_argument("--file_list2", "-2", default=None, type=str, metavar="/path/to/Condition2_rep1,/path/to/Condition2_rep2",
        help="Comma separated list of NanopolishComp files for label 2. (required if --sample_yaml not given)")
    parser_sc_sample_args.add_argument("--label1", type=str, metavar="Condition1", default="Condition1",
        help="Label for files in --file_list1 (default: %(default)s)")
    parser_sc_sample_args.add_argument("--label2", type=str, metavar="Condition2", default="Condition2",
        help="Label for files in --file_list2 (default: %(default)s)")
    parser_sc_io = parser_sc.add_argument_group('Input options')
    parser_sc_io.add_argument("--fasta", "-f", type=str, required=True,
        help="Fasta file used for mapping (required)")
    parser_sc_io.add_argument("--bed", type=str, default=None,
        help="BED file with annotation of transcriptome used for mapping (optional)")
    parser_sc_filtering = parser_sc.add_argument_group('Transcript filtering options')
    parser_sc_filtering.add_argument("--max_invalid_kmers_freq", type=float, default=0.1,
        help="Max fequency of invalid kmers (default: %(default)s)")
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

    # simreads subparser
    parser_sr = subparsers.add_parser('simreads', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
            Simulate reads as a NanopolishComp like file from a fasta file and an inbuild model\n
            * Minimal example without model alteration
                nanocompore simreads -f ref.fa -o results -n 50\n
            * Minimal example with alteration of model intensity loc parameter for 50% of the reads
            nanocompore simreads -f ref.fa -o results -n 50 --intensity_mod 2 --mod_reads_freq 0.5 --mod_bases_freq 0.2"""))
    parser_sr.set_defaults(func=simreads_main)
    parser_sr_io = parser_sr.add_argument_group('Input options')
    parser_sr_io.add_argument("--fasta", "-f", type=str, required=True,
        help="Fasta file containing references to use to generate artificial reads")
    parser_sr_modify = parser_sr.add_argument_group('Signal modification options')
    parser_sr_modify.add_argument("--intensity_mod", type=float, default=0,
        help="Fraction of intensity distribution SD by which to modify the intensity distribution loc value (default: %(default)s)")
    parser_sr_modify.add_argument("--dwell_mod", type=float, default=0,
        help="Fraction of dwell time distribution SD by which to modify the intensity distribution loc value (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_reads_freq", type=float, default=0,
        help="Frequency of reads to modify (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_freq", type=float, default=0.25,
        help="Frequency of bases to modify in each read (if possible) (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_type", type=str, default="A", choices=["A","T","C","G"],
        help="Base for which to modify the signal (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_extend_context", type=int, default=2,
        help="number of adjacent base affected by the signal modification following an harmonic series (default: %(default)s)")
    parser_sr_modify.add_argument("--min_mod_dist", type=int, default=6,
        help="Minimal distance between 2 bases to modify (default: %(default)s)")
    parser_sr_misc = parser_sr.add_argument_group('Other options')
    parser_sr_misc.add_argument("--run_type", type=str, default="RNA", choices=["RNA", "DNA"],
        help="Define the run type model to import (default: %(default)s)")
    parser_sr_misc.add_argument("--nreads_per_ref", "-n", type=int, default=100,
        help="Number of reads to generate per references (default: %(default)s)")
    parser_sr_misc.add_argument("--pos_rand_seed", type=int, default=42 ,
        help="Define a seed for randon position picking to get a deterministic behaviour (default: %(default)s)")
    parser_sr_misc.add_argument("--not_bound", action='store_true', default=False,
        help="Do not bind the values generated by the distributions to the observed min and max observed values from the model file (default: %(default)s)")

    # Eventalign_collapse subparser
    parser_ec = subparsers.add_parser("eventalign_collapse", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        Collapse the nanopolish eventalign output at kmers level and compute kmer level statistics
        * Minimal example
            nanocompore eventalign_collapse -i nanopolish_eventalign.tsv -outprefix out\n"""))
    parser_ec.set_defaults(func=eventalign_collapse_main)
    parser_ec_io = parser_ec.add_argument_group("Input options")
    parser_ec_io.add_argument("--eventalign", "-i", default=0,
        help="Path to a nanopolish eventalign tsv output file, or a list of file, or a regex (can be gzipped). It can be ommited if piped to standard input (default: piped to stdin)")
    parser_ec_rp = parser_ec.add_argument_group("Run parameters options")
    parser_ec_rp.add_argument("--n_lines", "-l", default=None , type=int ,
        help = "Number of lines to parse.(default: no limits")
    parser_ec_misc = parser_ec.add_argument_group("Other options")
    parser_ec_misc.add_argument("--nthreads", "-t", default=3, type=int,
        help="Total number of threads. 2 threads are reserved for the reader and the writer (default: %(default)s)")

    # Add common options for all parsers
    for sp in [parser_sc, parser_sr, parser_ec]:
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

    # Load eventalign_fn_dict from a YAML file or assemble eventalign_fn_dict for the command line option
    if args.sample_yaml:
        eventalign_fn_dict = args.sample_yaml
    elif args.file_list1 and args.file_list2:
        eventalign_fn_dict = build_eventalign_fn_dict(args.file_list1, args.file_list2, args.label1, args.label2)
    else:
        raise NanocomporeError("Samples eventalign files have to be provided with either `--sample_yaml` or `--file_list1` and `--file_list2`")

    # Init SampComp
    s = SampComp(
        eventalign_fn_dict = eventalign_fn_dict,
        max_invalid_kmers_freq = args.max_invalid_kmers_freq,
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
    db = s()

    # Save all reports
    if(db):
        db.save_all(pvalue_thr=args.pvalue_thr)

def simreads_main(args):
    """"""
    logger.warning("Running SimReads")

    # Run SimReads
    SimReads(
        fasta_fn = args.fasta,
        outpath = args.outpath,
        outprefix = args.outprefix,
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

def eventalign_collapse_main (args):
    """"""
    logger.warning("Running Eventalign_collapse")

    # Init Eventalign_collapse
    e = Eventalign_collapse (
        eventalign_fn = args.eventalign,
        outpath = args.outpath,
        outprefix = args.outprefix,
        overwrite = args.overwrite,
        n_lines = args.n_lines,
        nthreads = args.nthreads,
        progress = args.progress)

    # Run eventalign_collapse
    e()

#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    # execute only if run as a script
    main()
