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
from nanocompore.SimReads import SimReads
from nanocompore.common import NanocomporeError

# Main entry point
def main(args=None):
    parser = argparse.ArgumentParser(
            description="""
            Find signal differences between nanopolish eventalign collapsed files
            """)

    parser.add_argument('--version', '-v', action='version', version='v'+package_version)
    subparsers = parser.add_subparsers(help='Nanocompore implements the following subcommands', dest='sub-command')
    subparsers.required = True

    # Sampcomp subparser
    parser_sc = subparsers.add_parser('sampcomp', description="Compare 2 samples and find significant signal ")
    parser_sc.set_defaults(func=sampcomp_main)
    parser_sc_input = parser_sc.add_argument_group('Input options')
    parser_sc_input.add_argument("--file_list1", type=str, required=True, metavar="/path/to/Condition1_rep1,/path/to/Codition1_rep2", help="Comma separated list of NanopolishComp files for label 1 (required)")
    parser_sc_input.add_argument("--file_list2", type=str, required=True, metavar="/path/to/Condition2_rep1,/path/to/Codition2_rep2", help="Comma separated list of NanopolishComp files for label 2 (required)")
    parser_sc_input.add_argument("--label1", type=str, required=True, metavar="Condition1", default="Condition1", help="Label for files in --file_list1 (default: %(default)s)")
    parser_sc_input.add_argument("--label2", type=str, required=True, metavar="Condition2", default="Condition2", help="Label for files in --file_list2 (default: %(default)s)")
    parser_sc_input.add_argument("--fasta", "-f", type=str, required=True, help="Fasta file used for mapping (required)")
    parser_sc_input.add_argument("--bed", type=str, default=None, help="BED file with annotation of transcriptome used for mapping (optional)")
    parser_sc_output = parser_sc.add_argument_group('Output options')
    parser_sc_output.add_argument("--outpath", "-o", type=str, required=True, help="Path to the output folder (required)")
    parser_sc_output.add_argument("--overwrite", action='store_true', default=False, help="Use --outpath even if it exists already (default: %(default)s)")
    parser_sc_filtering = parser_sc.add_argument_group('Transcript filtering options')
    parser_sc_filtering.add_argument("--max_invalid_kmers_freq", type=float, default=0.1, help="Max fequency of invalid kmers (default: %(default)s)")
    parser_sc_filtering.add_argument("--min_coverage", type=int, default=50, help="Minimum coverage required in each condition to do the comparison (default: %(default)s)")
    parser_sc_filtering.add_argument("--downsample_high_coverage", type=int, default=None, help="Used for debug: transcripts with high covergage will be downsampled (default: %(default)s)")
    parser_sc_testing = parser_sc.add_argument_group('Statistical testing options')
    parser_sc_testing.add_argument("--comparison_methods", type=str, default="GMM,KS", help="Comma separated list of comparison methods. Valid methods are: GMM,KS,TT,MW. (default: %(default)s)")
    parser_sc_testing.add_argument("--sequence_context", type=int, default=2, choices=range(0,5), help="Sequence context for combining p-values (default: %(default)s)")
    parser_sc_testing.add_argument("--sequence_context_weights", type=str, default="uniform", choices=["uniform", "harmonic"], help="Type of weights to use for combining p-values")
    parser_sc_testing.add_argument("--pvalue_thr", type=float, default=0.05, help="Adjusted p-value threshold for reporting significant sites (default: %(default)s)")
    parser_sc_testing.add_argument("--logit", action='store_true', help="Use logistic regression testing also when all conditions have replicates (default: %(default)s)")
    parser_sc_testing.add_argument("--strict", type=bool, default=True, help="If True runtime warnings during the tests raise an error (default: %(default)s)")
    parser_sc_common = parser_sc.add_argument_group('Other options')
    parser_sc_common.add_argument("--nthreads", "-n", type=int, default=3, help="Number of threads (default: %(default)s)")
    parser_sc_common.add_argument("--log_level", type=str, default="info", choices=["warning", "info", "debug"], help="log level (default: %(default)s)")

    # Simulate_reads subparser
    parser_sr = subparsers.add_parser('simreads', description="Simulate reads in a NanopolishComp like file from a fasta file and an inbuild model")
    parser_sr.set_defaults(func=simreads_main)
    parser_sr_input = parser_sr.add_argument_group('Input options')
    parser_sr_input.add_argument("--fasta", "-f", type=str, required=True, help="Fasta file containing references to use to generate artificial reads (required)")
    parser_sr_input.add_argument("--run_type", type=str, default="RNA", help="Define the run type model to import (RNA or DNA) (default: %(default)s)")
    parser_sr_output = parser_sr.add_argument_group('Output options')
    parser_sr_output.add_argument("--outpath", "-o", type=str, default="./", help="Path to the output folder (default: %(default)s)")
    parser_sr_output.add_argument("--prefix", "-p", type=str, default="reads", help="text prefix for all the files generated by the function (default: %(default)s)")
    parser_sr_output.add_argument("--nreads_per_ref", type=int, default=100, help="Number of reads to generate per references (default: %(default)s)")
    parser_sr_modify = parser_sr.add_argument_group('Signal modification options')
    parser_sr_modify.add_argument("--intensity_mod_loc", type=float, default=0, help="value by which to modify the intensity distribution loc value (mode) (default: %(default)s)")
    parser_sr_modify.add_argument("--intensity_mod_scale", type=float, default=0 , help="value by which to modify the intensity distribution scale value (dispersion) (default: %(default)s)")
    parser_sr_modify.add_argument("--dwell_mod_loc", type=float, default=0, help="value by which to modify the dwell time distribution loc value (mode) (default: %(default)s)")
    parser_sr_modify.add_argument("--dwell_mod_scale", type=float, default=0, help="value by which to modify the dwell time distribution scale value (mode) (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_reads_freq", type=float, default=0, help="Frequency of reads to modify (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_freq", type=float, default=0, help="Frequency of bases to modify in each read (if possible) (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_bases_type", type=str, default="A", help="Base for which to modify the signal (default: %(default)s)")
    parser_sr_modify.add_argument("--mod_extend_context", type=int, default=2, help="number of adjacent base affected by the signal modification following an harmonic serries (default: %(default)s)")
    parser_sr_modify.add_argument("--min_mod_dist", type=int, default=6, help="Minimal distance between to bases to modify (default: %(default)s)")
    parser_sr_common = parser_sr.add_argument_group('Other options')
    parser_sr_common.add_argument("--pos_rand_seed", type=int, default=42 , help="Define a seed for randon position picking to get a deterministic behaviour (default: %(default)s)")
    parser_sr_common.add_argument("--distr_rand_seed", type=int, default=42 , help="Define a seed for randon distribution sampling to get a deterministic behaviour (default: %(default)s)")
    parser_sr_common.add_argument("--log_level", type=str, default="info", choices=["warning", "info", "debug"], help="log level (default: %(default)s)")

    # Downstream plot subparser
    parser_plot = subparsers.add_parser('plot', help="Run downstream analysis and plot results")
    parser_plot.set_defaults(func=plot)
    parser_plot.add_argument("--sampComp_db", type=str, help="path to SampCompDB")

    # Parse agrs and call subfunction
    args = parser.parse_args()
    args.func(args)

def sampcomp_main(args):

    # Check if output folder already exists
    outpath=Path(args.outpath)
    if outpath.exists():
        if not args.overwrite:
            raise NanocomporeError("%s already exists and --overwrite not specified" % args.outpath)
        elif not outpath.is_dir():
            raise NanocomporeError("%s is not a folder" % args.outpath )
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
        raise NanocomporeError("%s is not a valid file" % args.fasta )

    # Check if BED file exists
    if args.bed:
        bed_fn=Path(args.bed)
        if not bed_fn.is_file():
            raise NanocomporeError("%s is not a valid file" % args.bed )

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
        logit = args.logit,
        strict = args.strict,
        sequence_context = args.sequence_context,
        sequence_context_weights = args.sequence_context_weights,
        log_level = args.log_level)

    sc_out = s()

    #Save main report
    sc_out.save_report(output_fn="%s/nanocompore_results.txt" % outpath)
    sc_out.save_shift_stats(output_fn="%s/nanocompore_shift_stats.txt" % outpath)

    # Save bed and bedgraph files for each method used
    if args.bed:
        r = re.compile(".*_pvalue.*")
        methods = list(filter(r.match, list(sc_out.results)))
        out_bedpath = outpath / "bed_files"
        out_bedpath.mkdir(exist_ok=True)
        out_bedgpath = outpath / "bedgraph_files"
        out_bedgpath.mkdir(exist_ok=True)
        for m in methods:
            sc_out.save_to_bed(output_fn="%s/sig_sites_%s_thr%s.bed" %(out_bedpath, m, args.pvalue_thr), bedgraph=False, pvalue_field=m, pvalue_thr=args.pvalue_thr, span=5, title="Nanocompore Significant Sites")
            sc_out.save_to_bed(output_fn="%s/sig_sites_%s_thr%s.bedg" %(out_bedgpath, m, args.pvalue_thr), bedgraph=True, pvalue_field=m, title="Nanocompore Significant Sites")

def simreads_main(args):

    # Check if fasta file exists
    fasta_fn=Path(args.fasta)
    if not fasta_fn.is_file():
        raise NanocomporeError("%s is not a valid file" % args.fasta)

    SimReads (
        fasta_fn=args.fasta,
        outpath=args.outpath,
        prefix=args.prefix,
        run_type=args.run_type,
        nreads_per_ref=args.nreads_per_ref,
        intensity_mod_loc=args.intensity_mod_loc,
        intensity_mod_scale=args.intensity_mod_scale,
        dwell_mod_loc=args.dwell_mod_loc,
        dwell_mod_scale=args.dwell_mod_scale,
        mod_reads_freq=args.mod_reads_freq,
        mod_bases_freq=args.mod_bases_freq,
        mod_bases_type=args.mod_bases_type,
        mod_extend_context=args.mod_extend_context,
        min_mod_dist=args.min_mod_dist,
        pos_rand_seed=args.pos_rand_seed,
        distr_rand_seed=args.distr_rand_seed,
        log_level=args.log_level)

def plot(args):
    raise NanocomporeError("The plotting CLI methods haven't been implemented yet. Please load the the SampCompDB in jupyter for downstream analysis.")

def build_eventalign_fn_dict(labels, files):
    """ Build the eventalign_fn_dict
        labels: tuple of size 2 with sample lables
        files:  tuple of size 2 with comma separated lists of files
    """
    eventalign_fn_dict = dict()
    for cond in range(2):
        eventalign_fn_dict[ labels[cond] ] = { "%s%s" % (labels[cond], i): v for i,v in enumerate(files[cond].split(","), 1) }
    return(eventalign_fn_dict)


if __name__ == "__main__":
    # execute only if run as a script
    main()
