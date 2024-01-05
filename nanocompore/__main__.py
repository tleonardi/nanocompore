#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import argparse
import os
import shutil
import textwrap
import yaml

# Third party
from loguru import logger

# Local imports
from nanocompore import __version__ as package_version
from nanocompore import __description__ as package_description
from nanocompore.SampComp import SampComp
from nanocompore.Config import Config
from nanocompore.common import *

#~~~~~~~~~~~~~~MAIN PARSER ENTRY POINT~~~~~~~~~~~~~~#

def main(args=None):
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v'+package_version)

    subparser_description = textwrap.dedent("""
            nanocompore implements the following subcommands
            \t* template : Initialize a new input configuration file using the default template.
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
    logger.warning("Running SampComp")

    # Read the input config file
    with open(args.config, 'r') as f:
        config_file = yaml.safe_load(f)

    config = Config(config_file)

    # Check if output folder already exists
    try:
        mkdir(fn=config.get_outpath(), exist_ok=config.get_overwrite())
    except (NanocomporeError, FileExistsError) as E:
        raise NanocomporeError("Could not create the output folder. Try using `overwrite: true` in the input configuration or use another directory")


    # Init SampComp
    s = SampComp(config)

    # Run SampComp
    s()


def init_subcommand(args):
    """
    Initializes a new input configuration file using the default template.
    """
    template = os.path.join(os.path.dirname(__file__), 'template_config.yaml')
    if os.path.isfile(args.path) and not args.overwrite:
        raise NanocomporeError("Output file already exists. Use --overwrite to overwrite it.")
    shutil.copyfile(template, args.path)


#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    # execute only if run as a script
    main()
