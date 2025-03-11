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
from nanocompore.preprocessing import Uncalled4Preprocessor
from nanocompore.preprocessing import EventalignPreprocessor
from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.config import Config
from nanocompore.common import *


def main(args=None):
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v'+package_version)

    subparser_description = textwrap.dedent("""
            nanocompore implements the following subcommands
            \t* template : Initialize a new input configuration file using the default template.
            \t* preprocess : Preprocess the resquiggling data to prepare it for the subsequent analysis step.\n
            \t* eventalign_collapse : Parse eventalign data to process and store it to an intermediary efficient database for later analysis.\n
            \t* run : Compare 2 samples and find significant signal differences.\n""")
    subparsers = parser.add_subparsers(dest='subcommand',
                                       required=True,
                                       description=subparser_description)

    # preprocess subparser
    parser_sc = subparsers.add_parser('preprocess',
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=textwrap.dedent("Preprocess the resquiggling data to prepare it for comparison."))
    parser_sc.add_argument('config', type=str)
    parser_sc.set_defaults(func=preprocess_subcommand)

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

def preprocess_subcommand(args):
    """
    Runs the preprocessing step that will
    take the data from the chosen resquiggler
    and prepare it for the actual analysis.
    """
    logger.warning("Running the preprocessing step")

    # Read the input config file
    with open(args.config, 'r') as f:
        config_file = yaml.safe_load(f)

    try:
        config = Config(config_file)
    except Exception as e:
        msg = f"ERROR: {e}"
        sys.stderr.write(msg)
        exit(1)

    setup_logger(config, "preprocess.log")

    # Init the preprocessor
    if config.get_resquiggler() == "remora":
        RemoraPreprocessor(config)()
    elif config.get_resquiggler() == "uncalled4":
        preprocessor = Uncalled4Preprocessor(config)()
    elif config.get_resquiggler() == "eventalign":
        preprocessor = EventalignPreprocessor(config)()
    else:
        raise ArgumentError(f"Unsupported resquiggler {config.resquiggler}")


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

