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
    parser_sampComp = subparsers.add_parser('sampcomp', help="Run SumpComp")
    parser_sampComp.add_argument("--file_list1", type=str, help="Comma separated list of NanopolishComp files for label 1", required=True)
    parser_sampComp.add_argument("--file_list2", type=str, help="Comma separated list of NanopolishComp files for label 2", required=True)
    parser_sampComp.add_argument("--label1", type=str, help="Label for files in --file_list1", required=True)
    parser_sampComp.add_argument("--label2", type=str, help="Label for files in --file_list2", required=True)
    parser_sampComp.add_argument("--fasta", type=str, help="Fasta file used for mapping", required=True)
    parser_sampComp.add_argument("--bed", type=str, help="BED file with annotation of transcriptome used for mapping", default=None)
    parser_sampComp.add_argument("--comparison_methods", type=str, help="Comma separated list of comparison methods. Valid methods are: GMM,KS,TT,MW", default="GMM,KS")
    parser_sampComp.add_argument("--sequence_context", type=int, help="Sequence context for combining p-values", default=2, choices=range(0,5))
    parser_sampComp.add_argument("--outpath", "-o", type=str, help="Path to the output folder", required=True)
    parser_sampComp.add_argument("--overwrite", help="Use --outpath even if it exists already", action='store_true')
    parser_sampComp.add_argument("--max_invalid_kmers_freq", type=float, default=0.1, help="Max fequency of invalid kmers")
    parser_sampComp.add_argument("--min_coverage", type=int, default=50, help="Minimum coverage required in each condition to do the comparison")
    parser_sampComp.add_argument("--downsample_high_coverage", type=int, default=None, help="Used for debug: transcripts with high covergage will be downsampled")
    parser_sampComp.add_argument("--pvalue_thr", type=float, default=0.05, help="Adjusted p-value threshold for reporting significant sites")
    parser_sampComp.add_argument("--nthreads", "-n", type=int, default=3, help="Number of threads")
    parser_sampComp.add_argument("--loglevel", type=str, default="info", help="log level", choices=["warning", "info", "debug"])

    parser_sampComp.set_defaults(func=sample_compare_main)

    # Downstream subparser
    parser_downstream = subparsers.add_parser('downstream', help="Run downstream")
    parser_downstream.add_argument("--file", "-f", type=str, help="path to the summary file.")
    parser_downstream.set_defaults(func=sample_compare_main)

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
        logLevel = args.loglevel,
        sequence_context = args.sequence_context)

    sc_out = s()

    #Save main report
    sc_out.save_report(output_fn=f"{outpath}/nanocompore_results.txt")

    # Save bed and bedg files for each method used
    r = re.compile("adjusted_*")
    methods = list(filter(r.match, list(sc_out.results)))
    out_bedpath = outpath / "bed_files"
    out_bedpath.mkdir(exist_ok=True)
    out_bedgpath = outpath / "bedgraph_files"
    out_bedgpath.mkdir(exist_ok=True)
    for m in methods:
        sc_out.save_to_bed(output_fn=f"{out_bedpath}/sig_sites_{m}_thr{args.pvalue_thr}.bed", bedgraph=False, pvalue_field=m, pvalue_thr=args.pvalue_thr, span=5, title="Nanocompore Significant Sites")
        sc_out.save_to_bed(output_fn=f"{out_bedgpath}/sig_sites_{m}_thr{args.pvalue_thr}.bedg", bedgraph=True, pvalue_field=m, title="Nanocompore Significant Sites")

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
