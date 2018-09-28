# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, namedtuple
import shelve

# Third party
from pyfaidx import Fasta
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns

# Local package
from nanocompore.common import counter_to_str, access_file, NanocomporeError

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampCompDB (object):
    """ Wrapper over the shelve db"""

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#
    def __init__(self, db_fn, fasta_fn):
        # Check file

        for fn in (db_fn, fasta_fn):
            if not access_file (fn):
                raise NanocomporeError("Cannot access file {}".format(fn))

        # Try to get ref_id list from shelve db
        try:
            with shelve.open (db_fn, flag='r') as db:
                self.ref_id_list = list (db.keys())
                if not self.ref_id_list:
                    raise NanocomporeError("The result database is empty")
                self._db_fn = db_fn
        except:
            raise NanocomporeError("The result database cannot be opened")

        # Try to open Fasta file
        try:
            self._fasta = Fasta(fasta_fn)
        except:
            raise NanocomporeError("The fasta reference file cannot be opened")

    def __repr__ (self):
        """readable description of the object"""
        return "[{}] Number of references: {}\n".format(self.__class__.__name__, len(self))

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __len__ (self):
        return len (self.ref_id_list)

    def __iter__ (self):
        with shelve.open (self._db_fn, flag = "r") as db:
            for k, v in db.items():
                yield (k,v)

    def __getitem__(self, items):
        with shelve.open (self._db_fn, flag = "r") as db:
            if items in db:
                return db[items]
            else:
                return None

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot_pvalues (self, ref_id, start=None, end=None, pvalue_name="pvalue_median"):
        pass

    def plot_signal (self, ref_id, start, end):
        ref_fasta = self._fasta[ref_id]
        if start > end:
            raise NanocomporeError ("End coordinate has to be higher or equal to start")
        if start < 0:
            raise NanocomporeError ("Coordinates have to be higher that 0")
        if end > len(ref_fasta):
            raise NanocomporeError ("Coordinates have to be lower than the ref_id sequence length ({})".format(len(ref_fasta)))

        # Parse line position per position
        lt = namedtuple ("lt", ["pos", "sample", "median", "dwell"])
        l = []
        ref_pos_dict = self[ref_id]
        for pos in range (start, end+1):

            # Collect results for position
            if pos in ref_pos_dict:
                for sample in ("S1", "S2"):
                    for median, dwell in zip (ref_pos_dict[pos][sample+"_median"], ref_pos_dict[pos][sample+"_dwell"]):
                        l.append (lt (pos, sample, median, dwell))
            # If not coverage for position, just fill in with empty values
            else:
                for sample in ("S1", "S2"):
                    l.append (lt (pos, sample, None, None))

        # Create x label including the original sequence and its position
        x_lab = []
        for pos, base in zip (range (start, end+1), ref_fasta[start:end+1]):
            x_lab.append ("{}\n{}".format(pos, base))

        # Cast collected results to dataframe
        df = pd.DataFrame (l)
        fig, axes = pl.subplots(2, 1, figsize=(25,10))
        for variable, ax in zip (("median", "dwell"), axes):
            _  = sns.violinplot (x="pos", y=variable, hue="sample", data=df, split=True, ax=ax, inner="quartile", bw=0.75, linewidth=1)
            _ = ax.set_xticklabels(x_lab)
            _ = ax.set_xlabel ("")

        return (fig, axes)
