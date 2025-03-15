#!/usr/bin/env python3
import argparse
import os
import shutil
import textwrap
import yaml
import sys

from loguru import logger

from nanocompore import __version__ as package_version
from nanocompore import __description__ as package_description
from nanocompore.run import RunCmd
from nanocompore.preprocessing import RemoraPreprocessor
from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.config import Config
from nanocompore.common import Kit
from nanocompore.common import mkdir
from nanocompore.common import NanocomporeError


def main(args=None):
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v'+package_version)

    subparser_description = textwrap.dedent("""
            nanocompore implements the following subcommands
            \t* template : Initialize a new input configuration file using the default template.
            \t* eventalign_collapse : Parse eventalign data to process and store it to an intermediary efficient database for later analysis.\n
            \t* remora_resquiggle : Use Remora to resquiggle a sample and create an SQLite DB with the signal measurements.\n
            \t* run : Compare 2 samples and find significant signal differences.\n""")
    subparsers = parser.add_subparsers(dest='subcommand',
                                       required=True,
                                       description=subparser_description)

    # Sampcomp subparser
    parser_sc = subparsers.add_parser('run',
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=textwrap.dedent("Compare 2 samples and find significant signal differences."))
    parser_sc.add_argument('config', type=str)
    parser_sc.set_defaults(func=sampcomp_subcommand)

    # preprocess eventalign_collapse
    parser_sc = subparsers.add_parser('eventalign_collapse',
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=textwrap.dedent("Parse eventalign data to process and store it to an intermediary efficient database for later analysis."))
    parser_sc.add_argument('--ref', '-r', help="Transcriptome fasta reference.")
    parser_sc.add_argument('--input', '-i', help="Path to input eventalign file. If not provided, the input is read from stdin (useful for piping nanopolish/f5c eventalign directly).")
    parser_sc.add_argument('--output', '-o', help="Path to output SQLite database.")
    parser_sc.add_argument('--nthreads', '-n', help="Number of parallel processes to use for processing.", nargs='?', type=int, const=2, default=2)
    parser_sc.set_defaults(func=eventalign_collapse_subcommand)

    # preprocess remora_resquiggle
    parser_sc = subparsers.add_parser('remora_resquiggle',
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=textwrap.dedent("Use Remora to resquiggle a sample and create an SQLite DB with the signal measurements."))
    parser_sc.add_argument('--ref', '-r', help="Transcriptome fasta reference.", required=True)
    parser_sc.add_argument('--pod5', '-p', help="Path to input pod5 file containing the raw signal data.", required=True)
    parser_sc.add_argument('--bam', '-b', help="Path to input bam file containing the aligned reads.", required=True)
    parser_sc.add_argument('--output', '-o', help="Path to output SQLite database.", required=True)
    parser_sc.add_argument('--kit', '-k', help="Sequencing kit that was use (should be RNA002 or RNA004).", required=True)
    parser_sc.add_argument('--max-reads', '-m', help="Maximum number of reads to resquiggle for a transcript (default: 5000).", nargs='?', type=int, const=5000, default=5000)
    parser_sc.add_argument('--nthreads', '-n', help="Number of parallel processes to use for processing (default: 2).", nargs='?', type=int, const=2, default=2)
    parser_sc.set_defaults(func=remora_resquiggle_subcommand)

    # Init subparser
    parser_init = subparsers.add_parser('template',
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description=textwrap.dedent("Initialize a new input configuration file using the default template."))
    parser_init.add_argument('--overwrite', '-o', action='store_true', help="Overwrite existing config file.")
    parser_init.add_argument('path', type=str)
    parser_init.set_defaults(func=init_subcommand)

    # Parse agrs and
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # Call relevant subcommand function
    args.func(args)


#~~~~~~~~~~~~~~SUBCOMMAND FUNCTIONS~~~~~~~~~~~~~~#


def sampcomp_subcommand(args):
    """
    Runs the sample comparison subcommand.
    """

    # Read the input config file
    with open(args.config, 'r') as f:
        config_file = yaml.safe_load(f)

    try:
        config = Config(config_file)
    except Exception as e:
        msg = f"ERROR: {e}\n"
        sys.stderr.write(msg)
        exit(1)

    # If we'll show a progress bar then we
    # want to prevent the default loguru
    # logging on stderr.
    if config.get_progress:
        logger.remove()

    # Check if output folder already exists
    try:
        mkdir(fn=config.get_outpath(), exist_ok=config.get_result_exists_strategy() != 'stop')
    except (NanocomporeError, FileExistsError) as E:
        raise NanocomporeError("Could not create the output folder. Try using `overwrite: true` in the input configuration or use another directory")

    setup_logger(config, "run.log")

    run_cmd = RunCmd(config)
    run_cmd()


def eventalign_collapse_subcommand(args):
    """
    Parse eventalign data to process and store it
    to an intermediary efficient database for later
    analysis.
    """
    EventalignCollapser(args.input, args.ref, args.output, args.nthreads)()


def remora_resquiggle_subcommand(args):
    """
    Resquiggle a sample with Remora.
    """
    kit = Kit[args.kit]
    RemoraPreprocessor(args.ref,
                       args.pod5,
                       args.bam,
                       args.output,
                       kit,
                       args.max_reads,
                       args.nthreads)()


def init_subcommand(args):
    """
    Initializes a new input configuration file using the default template.
    """
    template = os.path.join(os.path.dirname(__file__), 'template_config.yaml')
    if os.path.isfile(args.path) and not args.overwrite:
        raise NanocomporeError("Output file already exists. Use --overwrite to overwrite it.")
    shutil.copyfile(template, args.path)


def setup_logger(config, file_name):
    if config.get_result_exists_strategy() == 'continue':
        logger_mode = 'a'
    else:
        logger_mode = 'w'
    logger.add(os.path.join(config.get_outpath(),
                            file_name),
                mode=logger_mode,
                level=config.get_log_level())


#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    main()

