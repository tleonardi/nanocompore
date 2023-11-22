#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Standard library imports
import argparse
from collections import *
import os
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
    parser.add_argument('config', type=str)

    # Parse agrs and
    args = parser.parse_args()

    # Read the input config file
    with open(args.config, 'r') as f:
        config_file = yaml.safe_load(f)

    config = Config(config_file)

    # Check if output folder already exists
    try:
        mkdir(fn=config.get_outpath(), exist_ok=config.get_overwrite())
    except (NanocomporeError, FileExistsError) as E:
        raise NanocomporeError("Could not create the output folder. Try using `overwrite: True` in the input configuration or use another directory")

    # Set logger
    log_fn = os.path.join(config.get_outpath(), f"{config.get_outprefix()}_nanocompore.log")

    set_logger(config.get_log_level(), log_fn=log_fn)

    # Run the comparison
    sampcomp_main(config)

#~~~~~~~~~~~~~~SUBCOMMAND FUNCTIONS~~~~~~~~~~~~~~#

def sampcomp_main(config):
    """"""
    logger.warning("Running SampComp")

    # Init SampComp
    s = SampComp(config)

    # Run SampComp
    s()
#~~~~~~~~~~~~~~CLI ENTRYPOINT~~~~~~~~~~~~~~#

if __name__ == "__main__":
    # execute only if run as a script
    main()
