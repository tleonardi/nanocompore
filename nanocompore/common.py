import datetime
import inspect
import multiprocessing as mp
import os
import sys
import time

from collections import OrderedDict
from collections import namedtuple
from enum import Enum

import numpy as np
import pysam

from jaxtyping import Float
from loguru import logger
from pyfaidx import Fasta


Kit = Enum('Kit', ['RNA002', 'RNA004'])

# Center positions taken from the Uncalled4 paper
# https://doi.org/10.1101/2024.03.05.583511 (Fig. 1a
# for RNA002 and Supplementary Fig. 4 for RNA004).
# The center position represents the base that
# has the largest influence on the measured current
# levels for the kmer. It's defined with a one-based
# index.
Kit.RNA002.len = 5
Kit.RNA002.center = 4
Kit.RNA004.len = 9
Kit.RNA004.center = 5


READ_ID_TYPE = np.uint32
SAMPLE_ID_TYPE = np.uint8
MEASUREMENTS_TYPE = np.float32


INTENSITY_POS = 0
DWELL_POS = 1
MOTOR_DWELL_POS = 2

EVENTALIGN = 'eventalign'
UNCALLED4 = 'uncalled4'
REMORA = 'remora'

TranscriptRow = namedtuple('TranscriptRow', 'ref_id id')


def get_pos_kmer(pos, seq, kit):
    """
    Takes a position <pos> (in 0-based transcriptomic coords)
    and the reference sequence and returns the k-mer that
    starts at the given position with k-mer length depending
    on the kit.
    """
    return seq[pos:pos+kit.len]


class NanocomporeError(Exception):
    """ Basic exception class for nanocompore module """
    def __init__(self, message=""):
        super().__init__(message)


class NanocomporeWarning(Warning):
    """ Basic Warning class for nanocompore module """
    pass


def log_init_state(loc):
    logger.debug("\ttimestamp: {}".format(str(datetime.datetime.now())))
    for i, j in loc.items():
        if type(j) in [int, float, complex, list, dict, str, bool, set, tuple]: # Avoid non standard types
            logger.debug("\t{}: {}".format(i,j))


def mkdir(fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs(fn, exist_ok=exist_ok)
    except Exception:
        raise NanocomporeError("Error creating output folder `{}`".format(fn))


def doc_func(func):
    """Parse the function description string"""

    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        if l.startswith("*"):
            break
        else:
            docstr_list.append(l)

    return "\n".join(docstr_list)


def make_arg_dict(func):
    """Parse the arguments default value, type and doc"""

    # Init method for classes
    if inspect.isclass(func):
        func = func.__init__

    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if p.name not in ["self","cls"]: # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if name not in ["kwargs","args"]: # Include but skip default required and type
                    # Get Annotation
                    if p.annotation != inspect._empty:
                        d[name]["type"] = p.annotation
                    # Get default value if available
                    if p.default == inspect._empty:
                        d[name]["required"] = True
                    else:
                        d[name]["default"] = p.default

        # Parse the docstring in a dict
        docstr_dict = OrderedDict()
        lab = None
        for l in inspect.getdoc(func).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    docstr_dict[lab] = []
                elif lab:
                    docstr_dict[lab].append(l)

        # Concatenate and copy doc in main dict
        for name in d.keys():
            if name in docstr_dict:
                d[name]["help"] = " ".join(docstr_dict[name])
        return d


def arg_opt(func, arg, **kwargs):
    """Get options corresponding to argumant name and deal with special cases"""
    arg_dict = make_arg_dict(func)[arg]

    if "default" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " (default: %(default)s)"

    if "type" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " [%(type)s]"

    # Special case for boolean args
    if arg_dict["type"] is bool:
        if arg_dict["default"] is False:
            arg_dict["action"] = 'store_true'
            del arg_dict["type"]
        elif arg_dict["default"] is True:
            arg_dict["action"] = 'store_false'
            del arg_dict["type"]

    # Special case for lists args
    elif arg_dict["type"] is list:
        arg_dict["nargs"] = '*'

    return arg_dict


def jhelp(f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    The docstring synthax should follow the same synthax as the one used for this function
    * f
        Function or method to display the help message for
    """
    # Private import as this is only needed if using jupyter
    from IPython.core.display import display, Markdown

    f_doc = doc_func(f)
    arg_doc = make_arg_dict(f)

    # Signature and function documentation
    s = "**{}** ({})\n\n{}\n\n---\n\n".format(f.__name__, ", ".join(arg_doc.keys()), f_doc)

    # Args doc
    for arg_name, arg_val in arg_doc.items():
        # Arg signature section
        s+= "* **{}**".format(arg_name)
        if "default" in arg_val:
            s+= " (default: {})".format(arg_val["default"])
        if "required" in arg_val:
            s+= " (required)"
        if "type" in arg_val:
            if type(list) == type:
                s+= " [{}]".format(arg_val["type"].__name__)
            else:
                s+= " [{}]".format(arg_val["type"])
        s+="\n\n"
        # Arg doc section
        if "help" in arg_val:
            s+= "{}\n\n".format(arg_val["help"])

    # Display in Jupyter
    display (Markdown(s))


def is_valid_fasta(file):
    """
    Returns a boolean indicating whether the given file is a valid FASTA file.
    """
    try:
        with Fasta(file) as fasta:
            # Fasta will return an empty iterator if the file is not properly
            # formatted FASTA file.
            return any(fasta)
    except IOError:
        # raise NanocomporeError("The fasta file cannot be opened")
        return False


def encode_kmer(kmer):
    encoding = 0
    for base in kmer:
        encoding <<= 2
        if base == "A" or base == "a":
            encoding |= 0
        elif base == "G" or base == "g":
            encoding |= 1
        elif base == "C" or base == "c":
            encoding |= 2
        elif base == "T" or base == "t":
            encoding |= 3
        else:
            raise RuntimeError(f"Bad nucleotide in kmer: {kmer}")
    return encoding


def decode_kmer(encoding, kmer_size):
    kmer = []
    for _ in range(kmer_size):
        base_code = encoding & ~(~0 << 2)
        base = ["A", "G", "C", "T"][base_code]
        kmer.append(base)
        encoding >>= 2
    return ''.join(kmer[::-1])


def monitor_workers(workers, delay_sec=5):
    """
    Check every *delay_sec* if any of the
    worker processes have been terminated
    with an error. If so, log an error and
    kill the all remaining workers and self.
    """
    terminated_ok = 0
    while terminated_ok < len(workers):
        time.sleep(delay_sec)
        for proc in workers:
            if proc.exitcode is not None:
                if proc.exitcode == 0:
                    terminated_ok += 1
                    continue
                logger.error(f"ERROR: A worker encountered an error (exitcode: {proc.exitcode}). "
                             "Will terminate all other workers and stop.")
                for child in mp.active_children():
                    child.terminate()
                sys.exit(1)


def get_references_from_bam(bam_path: str, has_data: bool=True) -> dict[str, int]:
    """
    Return references from the BAM's index.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.
    has_data : bool
        If set to true, return only references that have
        at least one mapped read to them.

    Returns
    -------
    dict[str, int]
        ref_id => number of reads
    """
    references = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for line in bam.get_index_statistics():
        if (has_data and line.mapped > 0) or not has_data:
            references[line.contig] = line.mapped
    return references


def get_reads_invalid_ratio(intensity: Float[np.ndarray, "positions reads"]): 
    """
    Calculate the ratio of missing kmers
    in the read. It takes the kmers with
    min and max position to determine the
    read length.

    This is the method employed for Uncalled4 and Remora.
    For eventalign there's a custom code (in eventalign collapse),
    that takes in consideration the richer information provided
    by the resquiggler.

    Parameters
    ----------
    intensity : npt.NDArray[np.float32]
        2D array with shape (positions, reads).
        The missing values should be np.nan. 
    Returns
    -------
    np.array
        Array with shape (reads,) containing
        the ratio of invalid positions for
        each read.
    """
    nreads = intensity.shape[1]
    ratios = np.empty(nreads)
    for n in range(nreads):
        valid_positions = np.where(~np.isnan(intensity[:, n]))[0]
        length = valid_positions.max() - valid_positions.min() + 1
        invalid_ratio = (length - len(valid_positions))/length
        ratios[n] = invalid_ratio
    return ratios

