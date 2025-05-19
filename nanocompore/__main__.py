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
from nanocompore.common import Kit
from nanocompore.common import NanocomporeError
from nanocompore.common import mkdir
from nanocompore.config import Config
from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.preprocessing import RemoraPreprocessor
from nanocompore import plotting
from nanocompore.run import RunCmd


def main(args=None):
    # General parser
    parser = argparse.ArgumentParser(description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', '-v', action='version', version='v'+package_version)

    subparser_description = textwrap.dedent(f"""
            nanocompore implements the following subcommands
            \t● {'template':<20} : Initialize a new input configuration file using the default template.
            \t● {'eventalign_collapse':<20} : Parse eventalign data to process and store it to an intermediary
                                           efficient SQLite database for later analysis.
            \t● {'remora_resquiggle':<20} : Use Remora to resquiggle a sample and create an intermediary
                                           efficient SQLite database for later analysis.
            \t● {'run':<20} : Compare 2 samples and find significant signal differences.
            \t● {'plot':<20} : Generate plots.""")
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
    parser_sc.add_argument('--ref', '-r', help="Transcriptome fasta reference.", required=True)
    parser_sc.add_argument('--input', '-i', help="Path to input eventalign file. If not provided, the input is read from stdin (useful for piping nanopolish/f5c eventalign directly).")
    parser_sc.add_argument('--output', '-o', help="Path to output SQLite database.", required=True)
    parser_sc.add_argument('--nthreads', '-n', help="Number of parallel processes to use for processing.", nargs='?', type=int, const=2, default=2)
    parser_sc.add_argument('--tmp', '-t', help="Directory where tmp files would be created (default: current directory).", nargs='?', type=str, const='.', default='.')
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

    parser_plot = subparsers.add_parser('plot',
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description=textwrap.dedent("Plotting functionality."))
    # plot_subparsers = parser_plot.add_subparsers(help='Plot type.')
    plot_subparsers_description = textwrap.dedent(f"""
            nanocompore plot implements the following subcommands
            \t● {'pvalues':<10} : Plot the p-values for a reference region.
            \t● {'signal':<10} : Plot the signal measurements for a region.
                                 Note: this will plot the data from the input
                                 files without applying any filtering performed
                                 by Nanocompore.
            \t● {'position':<10} : Plot the signal data for a specific position as
                                 a 2D plot. Note: this will plot the data as it is
                                 in the input files, without applying the filtering
                                 performed by Nanocompore.
            \t● {'gmm':<10} : Plot the GMM fitting for a position.
                                 This replicates all data filtering and produces
                                 identical GMM fitting to the one obtained by the
                                 run command.
            \t● {'coverage':<10} : Plot the read coverage over a reference region.
                                 Note: this will plot the data as it is in the input
                                 files, without applying the filtering performed
                                 by Nanocompore.""")
    plot_subparsers = parser_plot.add_subparsers(dest='subcommand',
                                                 required=True,
                                                 description=plot_subparsers_description)

    parser_pvalues = plot_subparsers.add_parser("pvalues",
                                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                                description=textwrap.dedent("Plot the p-values for a reference region."))
    parser_pvalues.add_argument('--config', '-c', type=str, required=True, help="Path to the input configuration YAML file.")
    parser_pvalues.add_argument('--output', '-o', type=str, required=True, help="Path and filename where the plot will be saved.")
    parser_pvalues.add_argument('reference', type=str, help="Reference name, matching a name from the FASTA reference used in the config.")
    parser_pvalues.add_argument('--start', type=int, help="0-based index on the reference.")
    parser_pvalues.add_argument('--end', type=int, help="0-based index on the reference.")
    parser_pvalues.add_argument('--threshold', type=float, help="Threshold p-value that will be drawn as a line.", default=0.01)
    parser_pvalues.add_argument('--kind', type=str, help="Kind of plot to draw. Default: lineplot", choices=['lineplot', 'barplot'], default='lineplot')
    parser_pvalues.add_argument('--figsize', type=str, help="Figure size. Default: 30,10", default='30,10')
    parser_pvalues.add_argument('--tests', type=str, help="Which tests to plot. By default all executed tests are shown.", default=None)
    parser_pvalues.add_argument('--palette', type=str, help="Color palette to use. Default: Dark2", default='Dark2')
    parser_pvalues.set_defaults(func=plot_pvalues)

    parser_signal = plot_subparsers.add_parser("signal",
                                               formatter_class=argparse.RawDescriptionHelpFormatter,
                                               description=textwrap.dedent("Plot the signal measurements for a region."))
    parser_signal.add_argument('--config', '-c', type=str, required=True, help="Path to the input configuration YAML file.")
    parser_signal.add_argument('--output', '-o', type=str, required=True, help="Path and filename where the plot will be saved.")
    parser_signal.add_argument('reference', type=str, help="Reference name, matching a name from the FASTA reference used in the config.")
    parser_signal.add_argument('--start', type=int, help="0-based index on the reference.")
    parser_signal.add_argument('--end', type=int, help="0-based index on the reference.")
    parser_signal.add_argument('--kind', type=str, help="Kind of plot to draw. Default: violinplot", choices=['violinplot', 'swarmplot', 'boxenplot'], default='violinplot')
    parser_signal.add_argument('--figsize', type=str, help="Figure size: WIDTH,HEIGHT. Default: 30,10", default='30,10')
    parser_signal.add_argument('--split_samples', action='store_true', help="Split results by sample instead of by condition.")
    parser_signal.add_argument('--markersize', type=int, help="Size of the points if swarmplot is used.", default=2)
    parser_signal.add_argument('--palette', type=str, help="Color palette to use. Default: Dark2", default='Dark2')
    parser_signal.set_defaults(func=plot_signal)

    parser_position = plot_subparsers.add_parser("position",
                                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                                 description=textwrap.dedent("Plot the signal data for a specific position as a 2D plot."))
    parser_position.add_argument('--config', '-c', type=str, required=True, help="Path to the input configuration YAML file.")
    parser_position.add_argument('--output', '-o', type=str, required=True, help="Path and filename where the plot will be saved.")
    parser_position.add_argument('reference', type=str, help="Reference name, matching a name from the FASTA reference used in the config.")
    parser_position.add_argument('position', type=int, help="Position to plot (0-based index).")
    parser_position.add_argument('--figsize', type=str, help="Figure size: WIDTH,HEIGHT. Default: 10,10", default='10,10')
    parser_position.add_argument('--point_size', type=int, help="Size of the data points", default=20)
    parser_position.add_argument('--xlim', type=str, help="Set specific range for the x-axis: MIN,MAX. By default it will be inferred from the data")
    parser_position.add_argument('--ylim', type=str, help="Set specific range for the x-axis: MIN,MAX. By default it will be inferred from the data")
    parser_position.add_argument('--kde', action='store_true', help="Plot the KDE of the intensity/dwell bivarariate distributions in the two samples.")
    parser_position.add_argument('--kde_levels', type=int, help="How many levels of the gaussian distributions to show", default=10)
    parser_position.add_argument('--palette', type=str, help="Use a single palette to select both point and gmm colors.", default='Dark2')
    parser_position.set_defaults(func=plot_position)

    parser_gmm = plot_subparsers.add_parser("gmm",
                                            formatter_class=argparse.RawDescriptionHelpFormatter,
                                            description=textwrap.dedent("Plot the GMM fitting for a single position."))
    parser_gmm.add_argument('--config', '-c', type=str, required=True, help="Path to the input configuration YAML file.")
    parser_gmm.add_argument('--output', '-o', type=str, required=True, help="Path and filename where the plot will be saved.")
    parser_gmm.add_argument('reference', type=str, help="Reference name, matching a name from the FASTA reference used in the config.")
    parser_gmm.add_argument('position', type=int, help="Position to plot (0-based index).")
    parser_gmm.add_argument('--figsize', type=str, help="Figure size: WIDTH,HEIGHT. Default: 10,10", default='10,10')
    parser_gmm.add_argument('--point_size', type=int, help="Size of the data points", default=20)
    parser_gmm.add_argument('--xlim', type=str, help="Set specific range for the x-axis: MIN,MAX. By default it will be inferred from the data")
    parser_gmm.add_argument('--ylim', type=str, help="Set specific range for the x-axis: MIN,MAX. By default it will be inferred from the data")
    parser_gmm.add_argument('--gmm_levels', type=int, help="How many levels of the gaussian distributions to show", default=4)
    parser_gmm.add_argument('--palette', type=str, help="Use a single palette to select both point and gmm colors.", default='Dark2')
    parser_gmm.add_argument('--point_palette', type=str, help="Which palette to use for plotting the points.")
    parser_gmm.add_argument('--gmm_palette', type=str, help="Which palette to use for plotting the gaussian distributions.")
    parser_gmm.set_defaults(func=plot_gmm)

    parser_coverage = plot_subparsers.add_parser("coverage",
                                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                                 description=textwrap.dedent("Plot the read coverage for a region."))
    parser_coverage.add_argument('--config', '-c', type=str, required=True, help="Path to the input configuration YAML file.")
    parser_coverage.add_argument('--output', '-o', type=str, required=True, help="Path and filename where the plot will be saved.")
    parser_coverage.add_argument('reference', type=str, help="Reference name, matching a name from the FASTA reference used in the config.")
    parser_coverage.add_argument('--start', type=int, help="0-based index on the reference.")
    parser_coverage.add_argument('--end', type=int, help="0-based index on the reference.")
    parser_coverage.add_argument('--figsize', type=str, help="Figure size: WIDTH,HEIGHT. Default: 30,10", default='12,4')
    parser_coverage.add_argument('--split_samples', action='store_true', help="Split results by sample instead of by condition.")
    parser_coverage.add_argument('--palette', type=str, help="Use a single palette to select both point and gmm colors.", default='Dark2')
    parser_coverage.set_defaults(func=plot_coverage)

    if len(sys.argv) == 1:
        # no argument passed
        args = parser.parse_args(['--help'])
    elif len(sys.argv) == 2 and sys.argv[1] == 'plot':
        # plot called without arguments
        args = parser.parse_args(['plot', '--help'])
    else:
        args = parser.parse_args()

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
    EventalignCollapser(args.input, args.ref, args.output, args.nthreads, args.tmp)()


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


def plot_pvalues(args: argparse.Namespace):
    """
    Plot pvalues along a reference region.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
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

    figsize = tuple(map(int, args.figsize.split(',')))
    tests = args.tests.split(',') if args.tests else None

    fig = plotting.plot_pvalues(config,
                                args.reference,
                                args.start,
                                args.end,
                                args.kind,
                                args.threshold,
                                figsize,
                                tests,
                                args.palette)
    fig.savefig(args.output)


def plot_signal(args: argparse.Namespace):
    """
    Plot the signal measurements along a reference region.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
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

    figsize = tuple(map(int, args.figsize.split(',')))

    fig = plotting.plot_signal(config,
                               args.reference,
                               args.start,
                               args.end,
                               args.kind,
                               figsize,
                               args.split_samples,
                               args.markersize,
                               args.palette)
    fig.savefig(args.output)


def plot_position(args: argparse.Namespace):
    """
    Plot the signal measurements for a specific position.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
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

    figsize = tuple(map(int, args.figsize.split(',')))
    xlim = None
    if args.xlim is not None:
        xlim = tuple(map(int, args.xlim.split(',')))
    ylim = None
    if args.ylim is not None:
        ylim = tuple(map(int, args.ylim.split(',')))

    fig = plotting.plot_position(config,
                                 args.reference,
                                 args.position,
                                 figsize,
                                 args.point_size,
                                 xlim,
                                 ylim,
                                 args.kde,
                                 args.kde_levels,
                                 args.palette)
    fig.savefig(args.output)

def plot_gmm(args: argparse.Namespace):
    """
    Plot GMM for a position.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
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

    figsize = tuple(map(int, args.figsize.split(',')))
    xlim = None
    if args.xlim is not None:
        xlim = tuple(map(int, args.xlim.split(',')))
    ylim = None
    if args.ylim is not None:
        ylim = tuple(map(int, args.ylim.split(',')))

    fig = plotting.plot_gmm(config,
                            args.reference,
                            args.position,
                            figsize,
                            args.point_size,
                            xlim,
                            ylim,
                            args.gmm_levels,
                            args.palette,
                            args.point_palette,
                            args.gmm_palette)
    fig.savefig(args.output)


def plot_coverage(args: argparse.Namespace):
    """
    Plot the read coverage for a region.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
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

    figsize = tuple(map(int, args.figsize.split(',')))

    fig = plotting.plot_coverage(config,
                                 args.reference,
                                 args.start,
                                 args.end,
                                 figsize,
                                 args.split_samples,
                                 args.palette)
    fig.savefig(args.output)


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

