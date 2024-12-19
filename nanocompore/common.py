# -*- coding: utf-8 -*-

from email.policy import default
from enum import Enum
import sys
import os
from collections import *
import inspect
import datetime

import numpy as np

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
Kit.RNA004.center = 6


def is_valid_position(pos, seq_len, kit):
    """
    Takes a position (in 0-based transcriptomic coords),
    the length of the reference transcript and the kit
    used and returns whether the position is within the
    trimmed range of the transcript.
    """
    # E.g. if we have 5mer model with center at base 4 and
    # a transcript of length 100 then positions 0, 1, and 2
    # on the left and 99 on the right would be invalid.
    if pos < kit.center - 1 or pos >= seq_len - (kit.len - kit.center):
        return False
    return True


def get_pos_kmer(pos, seq, kit):
    """
    Takes a position <pos> (in 0-based transcriptomic coords)
    and the reference sequence and returns the kmer that
    has its central (most influential) base at <pos>, depending
    on the kit.
    E.g. for RNA002, where the 4th base of the 5mer is the center
    it would return the seq[pos-3:pos+2].
    """
    return seq[(pos - kit.center + 1):(pos + kit.len - kit.center + 1)]


class NanocomporeError (Exception):
    """ Basic exception class for nanocompore module """
    def __init__(self, message=""):
        logger.error(message)
        super().__init__(message)


class NanocomporeWarning (Warning):
    """ Basic Warning class for nanocompore module """
    pass


def build_eventalign_fn_dict(pod5_list1, bam_list1, pod5_list2, bam_list2, label1, label2):
    """
    Build the eventalign_fn_dict from file lists and labels
    """

    pod5_list1 = pod5_list1.split(",")
    bam_list1 = bam_list1.split(",")

    pod5_list2 = pod5_list2.split(",")
    bam_list2 = bam_list2.split(",")

    if len(pod5_list1) == len(bam_list1) and len(pod5_list2) == len(bam_list2):
        d = defaultdict(dict)
        d[label1] = build_condition_dict(pod5_list1, bam_list1, label1)
        d[label2] = build_condition_dict(pod5_list2, bam_list2, label2)
        return d

    else:
        pod5_list1_len = len(pod5_list1)
        bam_list1_len = len(bam_list1)
        pod5_list2_len = len(pod5_list2)
        bam_list2_len = len(bam_list2)
        raise NanocomporeError (f"Mismatch in file list lengths:\npod5_list1 {pod5_list1_len}; bam_list1 {bam_list1_len}\npod5_list2 {pod5_list2_len}; bam_list2 {bam_list2_len}\n")


def build_condition_dict(pod5_list, bam_list, label):
    condition_list = defaultdict()
    for i, (pod5, bam) in enumerate(zip(pod5_list, bam_list)):
        condition_list[f"{label}_{i}"] = {'pod5':pod5, 'bam':bam}
    return condition_list


def build_eventalign_fn_dict_from_tsv(infile):
    fn_dict = defaultdict(dict)
    num_samples = set()
    with open(infile, 'r') as tsv:
        for i, line in enumerate(tsv):
            if line:
                try:
                    sample, cond, pod5, bam = line.strip().split('\t')
                    num_samples.add(sample)
                    fn_dict[cond][sample] = {'pod5':pod5, 'bam':bam}
                except:
                    raise NanocomporeError (f"Error with entry {i} in input samples tsv file\n")

    if len(num_samples) != i+1:
        raise NanocomporeError (f"Not all sample labels are unique\nCheck sample label column in input samples tsv file\n")

    return fn_dict


def set_logger (log_level, log_fn=None):
    log_level = log_level.upper()
    logger.remove()
    logger.add(sys.stderr, format="{time} {level} - {process.name} | {message}", enqueue=True, level=log_level)
    if log_fn:
        if os.path.isfile(log_fn):
            os.remove(log_fn)
        logger.add(log_fn, format="{time} {level} - {process.name} | {message}", enqueue=True, level="DEBUG")


def log_init_state(loc):
    #logger.debug("\tpackage_name: {}".format(pkg.__name__))
    #logger.debug("\tpackage_version: {}".format(pkg.__version__))
    logger.debug("\ttimestamp: {}".format(str(datetime.datetime.now())))
    for i, j in loc.items():
        if type(j) in [int, float, complex, list, dict, str, bool, set, tuple]: # Avoid non standard types
            logger.debug("\t{}: {}".format(i,j))


def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanocomporeError ("Error creating output folder `{}`".format(fn))


def access_file (fn, **kwargs):
    """ Check if the file is readable """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)


def numeric_cast_list (l):
    """ Cast str values to integer or float from a list """
    l2 = []
    for v in l:
        l2.append(numeric_cast(v))
    return l2


def numeric_cast_dict (keys, values):
    """ Cast str values to integer or float from a list """
    d = OrderedDict()
    for k, v in zip(keys, values):
        d[k] = numeric_cast(v)
    return d


def numeric_cast (v):
    if type(v)== str:
        try:
            v = int(v)
        except ValueError:
            try:
                v = float(v)
            except ValueError:
                pass
    return v


def counter_to_str (c):
    """ Transform a counter dict to a tabulated str """
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m


def all_values_in (required_val_list, all_val_list):
    """Check that all values in required_val_list are found in all_val_list"""
    for v in required_val_list:
        if not v in all_val_list:
            return False
    return True


def doc_func (func):
    """Parse the function description string"""

    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        if l.startswith("*"):
            break
        else:
            docstr_list.append(l)

    return "\n".join(docstr_list)


def make_arg_dict (func):
    """Parse the arguments default value, type and doc"""

    # Init method for classes
    if inspect.isclass(func):
        func = func.__init__

    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if not p.name in ["self","cls"]: # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if not name in ["kwargs","args"]: # Include but skip default required and type
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
        lab=None
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


def arg_opt (func, arg, **kwargs):
    """Get options corresponding to argumant name and deal with special cases"""
    arg_dict = make_arg_dict(func)[arg]

    if "default" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " (default: %(default)s)"

    if "type" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " [%(type)s]"

    # Special case for boolean args
    if arg_dict["type"] == bool:
        if arg_dict["default"] == False:
            arg_dict["action"] = 'store_true'
            del arg_dict["type"]
        elif arg_dict["default"] == True:
            arg_dict["action"] = 'store_false'
            del arg_dict["type"]

    # Special case for lists args
    elif arg_dict["type"] == list:
        arg_dict["nargs"]='*'

    return arg_dict


def jhelp (f:"python function or method"):
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

DROP_METADATA_TABLE_QUERY = """
DROP TABLE IF EXISTS metadata;
"""
CREATE_METADATA_TABLE_QUERY = """
CREATE TABLE metadata (
    key TEXT NOT NULL,
    value TEXT NOT NULL,
    PRIMARY KEY(key)
);
"""
DROP_KMER_DATA_TABLE_QUERY = """
DROP TABLE IF EXISTS kmer_data;
"""
CREATE_KMER_DATA_TABLE_QUERY = """
CREATE TABLE kmer_data (
    transcript_id INTEGER NOT NULL,
    pos INTEGER NOT NULL,
    kmer TEXT NOT NULL,
    samples BLOB NOT NULL,
    reads BLOB NOT NULL,
    intensity BLOB,
    intensity_std BLOB,
    dwell BLOB
);
"""
CREATE_INTERMEDIARY_KMER_DATA_TABLE_QUERY = """
CREATE TABLE kmer_data (
    transcript_id INTEGER NOT NULL,
    pos INTEGER NOT NULL,
    kmer TEXT NOT NULL,
    reads BLOB NOT NULL,
    intensity BLOB,
    intensity_std BLOB,
    dwell BLOB
);
"""
INSERT_KMER_DATA_QUERY = """
INSERT INTO kmer_data VALUES (?, ?, ?, ?, ?, ?, ?, ?);
"""
INSERT_INTERMEDIARY_KMER_DATA_QUERY = """
INSERT INTO kmer_data VALUES (?, ?, ?, ?, ?, ?, ?);
"""
DROP_READS_TABLE_QUERY = """
DROP TABLE IF EXISTS reads;
"""
CREATE_READS_TABLE_QUERY = """
CREATE TABLE reads (
    read TEXT NOT NULL,
    id INTEGER NOT NULL,
    invalid_kmers REAL
);
"""
INSERT_READS_QUERY = """
INSERT INTO reads VALUES(?, ?, ?);
"""
DROP_READS_ID_INDEX_QUERY = """
DROP INDEX IF EXISTS reads_id_index;
"""
CREATE_READS_ID_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS reads_id_index ON reads (id);
"""
DROP_KMERS_INDEX_QUERY = """
DROP INDEX IF EXISTS kmer_index;
"""
CREATE_KMERS_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS kmer_index ON kmer_data (transcript_id, pos);
"""

DROP_TRANSCRIPTS_TABLE_QUERY = """
DROP TABLE IF EXISTS transcripts;
"""
CREATE_TRANSCRIPTS_TABLE_QUERY = """
CREATE TABLE transcripts (
    reference TEXT NOT NULL,
    id INTEGER NOT NULL,
    PRIMARY KEY(reference)
);
"""
DROP_TRANSCRIPTS_ID_INDEX_QUERY = """
DROP INDEX IF EXISTS transcripts_id_index;
"""
CREATE_TRANSCRIPTS_ID_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS transcripts_id_index ON transcripts (id);
"""
INSERT_TRANSCRIPTS_QUERY = """
INSERT INTO transcripts VALUES(?, ?);
"""

READ_ID_TYPE = np.uint32
SAMPLE_ID_TYPE = np.uint8
REMORA_MEASUREMENT_TYPE = np.float32
UNCALLED4_MEASUREMENT_TYPE = np.int16
EVENTALIGN_MEASUREMENT_TYPE = np.float32

DB_METADATA_RESQUIGGLER_KEY = 'resquiggler'
DB_METADATA_READ_ID_TYPE_KEY = 'read_id_type'
DB_METADATA_SAMPLE_ID_TYPE_KEY = 'sample_id_type'
DB_METADATA_MEASUREMENT_TYPE_KEY = 'measurement_type'
DB_METADATA_SAMPLE_LABELS_KEY = 'sample_labels'
DB_METADATA_CONDITION_SAMPLES_KEY = 'condition_samples'


def get_measurement_type(resquiggler):
    if resquiggler == 'remora':
        return REMORA_MEASUREMENT_TYPE
    elif resquiggler == 'uncalled4':
        return UNCALLED4_MEASUREMENT_TYPE
    elif resquiggler == 'eventalign':
        return EVENTALIGN_MEASUREMENT_TYPE
    else:
        raise NotImplementedError(f"Unsupported resquiggler '{self._config.get_resquiggler()}'")


class Indexer:
    def __init__(self, initial_index=1):
        self._ids = {}
        self._current_id = initial_index


    def add(self, elements):
        """
        Add elements to the indexer.

        All new elements will be assigned a unique index
        and the new mappings will be returned to the
        caller as a list of (element, id) pairs.
        """
        new_elements = [element
                        for element in elements
                        if element not in self._ids]
        new_mappings = []
        for element in new_elements:
            self._ids[element] = self._current_id
            new_mappings.append((element, self._current_id))
            self._current_id += 1
        return new_mappings


    def get_ids(self, elements):
        return [self._ids[element] for element in elements]


    @property
    def current_id(self):
        return self._current_id


TranscriptRow = namedtuple('TranscriptRow', 'ref_id id')

