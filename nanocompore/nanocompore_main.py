#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys
from pathlib import Path
from collections import OrderedDict
import re
import json


# Local imports
from nanocompore import __version__ as package_version
from nanocompore import __name__ as package_name
from nanocompore.SampComp import SampComp
from nanocompore.Whitelist import Whitelist
from nanocompore.common import NanocomporeError


def main(args=None):
    parser = argparse.ArgumentParser(
            description="""
            Find signal differences between nanopolish eventalign collapsed files 
            """)

    parser.add_argument('--version', '-v', action='version', version='v'+package_version)
    subparsers = parser.add_subparsers(help='Nanocompore implements the following subcommands', dest='sub-command')
    subparsers.required = True
   
    # Sampcomp subparser
    parser_sampComp = subparsers.add_parser('sampcomp', help="Run SumpComp", formatter_class=argparse.MetavarTypeHelpFormatter, add_help=False)
    parser_sampComp.set_defaults(func=sample_compare_main)

    parser_input = parser_sampComp.add_argument_group('Input files')
    parser_input.add_argument("--file_list1", type=str, required=True, metavar="/path/to/Condition1_rep1,/path/to/Codition1_rep2", help="Comma separated list of NanopolishComp files for label 1 (required)")
    parser_input.add_argument("--file_list2", type=str, required=True, metavar="/path/to/Condition2_rep1,/path/to/Codition2_rep2", help="Comma separated list of NanopolishComp files for label 2 (required)")
    parser_input.add_argument("--label1", type=str, required=True, metavar="Condition1", default="Condition1", help="Label for files in --file_list1 (default: %(default)s)")
    parser_input.add_argument("--label2", type=str, required=True, metavar="Condition2", default="Condition2", help="Label for files in --file_list2 (default: %(default)s)")
    parser_input.add_argument("--fasta", type=str, required=True, help="Fasta file used for mapping (required)")
    parser_input.add_argument("--bed", type=str, default=None, help="BED file with annotation of transcriptome used for mapping (optional)")

    parser_output = parser_sampComp.add_argument_group('Output paths')
    parser_output.add_argument("--outpath", "-o", type=str, required=True, help="Path to the output folder (required)")
    parser_output.add_argument("--overwrite", action='store_true', help="Use --outpath even if it exists already (default: %(default)s)")

    parser_filtering = parser_sampComp.add_argument_group('Transcript filtering')
    parser_filtering.add_argument("--max_invalid_kmers_freq", type=float, default=0.1, help="Max fequency of invalid kmers (default: %(default)s)")
    parser_filtering.add_argument("--min_coverage", type=int, default=50, help="Minimum coverage required in each condition to do the comparison (default: %(default)s)")
    parser_filtering.add_argument("--downsample_high_coverage", type=int, default=None, help="Used for debug: transcripts with high covergage will be downsampled (default: %(default)s)")

    parser_testing = parser_sampComp.add_argument_group('Statistical testing')
    parser_testing.add_argument("--comparison_methods", type=str, default="GMM,KS", help="Comma separated list of comparison methods. Valid methods are: GMM,KS,TT,MW. (default: %(default)s)")
    parser_testing.add_argument("--sequence_context", type=int, default=2, choices=range(0,5), help="Sequence context for combining p-values (default: %(default)s)")
    parser_testing.add_argument("--sequence_context_weights", type=str, default="uniform", choices=["uniform", "harmonic"], help="Type of weights to use for combining p-values")
    parser_testing.add_argument("--pvalue_thr", type=float, default=0.05, help="Adjusted p-value threshold for reporting significant sites (default: %(default)s)")
    parser_output.add_argument("--force_logit", action='store_true', help="Use logistic regression testing even if all conditions have replicates (default: %(default)s)")

    parser_common = parser_sampComp.add_argument_group('Other options')
    parser_common.add_argument("--nthreads", "-n", type=int, default=3, help="Number of threads (default: %(default)s)")
    parser_common.add_argument("--loglevel", type=str, default="info", choices=["warning", "info", "debug"], help="log level (default: %(default)s)")
    parser_common.add_argument("--help", "-h", action="help", help="Print this help message")



    # Downstream subparser
    parser_plot = subparsers.add_parser('plot', help="Run downstream")
    parser_plot.add_argument("--sampComp_db", type=str, help="path to SampCompDB")
    parser_plot.set_defaults(func=plot)

    args = parser.parse_args()
    args.func(args)

def sample_compare_main(args):

    # Check if output folder already exists
    outpath=Path(args.outpath)
    if outpath.exists():
        if not args.overwrite:
            raise NanocomporeError(f"{args.outpath} already exists and --overwrite not specified")
        elif not outpath.is_dir():
            raise NanocomporeError(f"{args.outpath} is not a folder")
    else:
        outpath.mkdir(parents=True, exist_ok=False)

    # Save command line arguments to file 
    sampcomp_log=outpath/'sample_compare.log'
    with sampcomp_log.open(mode='w') as f:
        json.dump({k:v for k,v in vars(args).items() if k!="func" }, f, indent=2)

    # Assemble eventalign_fn_dict
    eventalign_fn_dict = build_eventalign_fn_dict((args.label1, args.label2), (args.file_list1, args.file_list2))

    # Check if fasta file exists
    fasta_fn=Path(args.fasta)
    if not fasta_fn.is_file():
        raise NanocomporeError(f"{args.fasta} is not a valid file")

    # Check if BED file exists
    if args.bed:
        bed_fn=Path(args.bed)
        if not bed_fn.is_file():
            raise NanocomporeError(f"{args.bed} is not a valid file")

    # Check at least 3 threads
    if args.nthreads < 3:
        raise NanocomporeError("The minimum number of threads is 3")

    s = SampComp (
        eventalign_fn_dict = eventalign_fn_dict,
        max_invalid_kmers_freq=args.max_invalid_kmers_freq,
        output_db_fn = args.outpath+"/sampCompDB.db",
        fasta_fn = args.fasta,
        bed_fn = args.bed,
        nthreads = args.nthreads,
        min_coverage = args.min_coverage,
        downsample_high_coverage = args.downsample_high_coverage,
        comparison_method = args.comparison_methods,
        force_logit = args.force_logit,
        logLevel = args.loglevel,
        sequence_context = args.sequence_context,
        sequence_context_weights = args.sequence_context_weights)

    sc_out = s()

    #Save main report
    sc_out.save_report(output_fn=f"{outpath}/nanocompore_results.txt")
    sc_out.save_shift_stats(output_fn=f"{outpath}/nanocompore_shift_stats.txt")

    # Save bed and bedg files for each method used
    if args.bed:
        r = re.compile(".*_pvalue.*")
        methods = list(filter(r.match, list(sc_out.results)))
        out_bedpath = outpath / "bed_files"
        out_bedpath.mkdir(exist_ok=True)
        out_bedgpath = outpath / "bedgraph_files"
        out_bedgpath.mkdir(exist_ok=True)
        for m in methods:
            sc_out.save_to_bed(output_fn=f"{out_bedpath}/sig_sites_{m}_thr{args.pvalue_thr}.bed", bedgraph=False, pvalue_field=m, pvalue_thr=args.pvalue_thr, span=5, title="Nanocompore Significant Sites")
            sc_out.save_to_bed(output_fn=f"{out_bedgpath}/sig_sites_{m}_thr{args.pvalue_thr}.bedg", bedgraph=True, pvalue_field=m, title="Nanocompore Significant Sites")

def plot(args):
    raise NanocomporeError("The plotting CLI methods haven't been implemented yet. Please load the the SampCompDB in jupyter for downstream analysis.")

def build_eventalign_fn_dict(labels, files):
    """ Build the eventalign_fn_dict
        labels: tuple of size 2 with sample lables
        files:  tuple of size 2 with comma separated lists of files
    """
    eventalign_fn_dict = dict()
    for cond in range(2):
        eventalign_fn_dict[ labels[cond] ] = { f"{labels[cond]}{i}": v for i,v in enumerate(files[cond].split(","), 1) }
    return(eventalign_fn_dict)


if __name__ == "__main__":
    # execute only if run as a script
    main()
