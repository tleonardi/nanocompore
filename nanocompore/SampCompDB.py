# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, namedtuple
import shelve

# Third party
from pyfaidx import Fasta
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns
from bedparse import bedline
from statsmodels.stats.multitest import multipletests

# Local package
from nanocompore.common import counter_to_str, access_file, NanocomporeError

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampCompDB (object):
    """ Wrapper over the result shelve SampComp """

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#
    def __init__(self, db_fn, fasta_fn):
        """
        Import a shelve db and a fasta reference file. Automatically returned by SampComp
        Can also be manually created from an existing shelve db output
        db_fn: Path where to write the result database
        fasta_fn: Path to a fasta file corresponding to the reference used for read alignemnt
        """
        # Check file
        for fn in (db_fn, fasta_fn):
            if not access_file (fn):
                raise NanocomporeError("Cannot access file {}".format(fn))

        # Try to get ref_id list from shelve db
        try:
            with shelve.open (db_fn, flag='r') as db:
                self.ref_id_list = list (db.keys())
                try: 
                    metadata = db['__metadata']
                except KeyError:
                    raise NanocomporeError("The result database does not contain metadata")
                self._comparison_method=metadata['comparison_method']
                self._sequence_context=metadata['sequence_context']
                self.ref_id_list.remove('__metadata')
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

    ################################### TO DO ##################################
    def to_bed (self, output_fn):
        pass

    def results(self, bed_fn=None, adjust=True):
        # Compose a lists with the name of the results
        tests=[]
        for method in self._comparison_method:
            if method in ["mann_whitney", "MW"]:
                tests.append("pvalue_mann_whitney_median")
                tests.append("pvalue_mann_whitney_dwell")
            elif method in ["kolmogorov_smirnov", "KS"]:
                tests.append("pvalue_kolmogorov_smirnov_median")
                tests.append("pvalue_kolmogorov_smirnov_dwell")
            elif method in ["t_test", "TT"]:
                tests.append("pvalue_t_test_median")
                tests.append("pvalue_t_test_dwell")
            elif method in ["kmean"]:
                tests.append("pvalue_kmeans")
        if self._sequence_context:
            c=str(self._sequence_context)
            tests+=[t+"_context="+c for t in tests]

        df = pd.DataFrame([dict({a:b for a,b in v.items() if a in tests}, ref=ref_id, pos=k)  for ref_id in self.ref_id_list for k,v in self[ref_id].items() ]).fillna(1)

        if bed_fn:
            bed_annot={}
            try:
                with open(bed_fn) as tsvfile:
                    for line in tsvfile:
                        record_name=line.split('\t')[3]
                        if( record_name in self.ref_id_list):
                            bed_annot[record_name]=bedline(line.split('\t'))
                if len(bed_annot) != len(self.ref_id_list):
                    raise NanocomporeError("Some references are missing from the BED file provided")
            except:
                raise NanocomporeError("Can't open BED file")

            df['genomicPos'] = df.apply(lambda row: bed_annot[row['ref']].tx2genome(coord=row['pos']),axis=1)
            df=df[['ref', 'pos', 'genomicPos']+tests]
        if adjust:
            for col in tests:
                df['adjusted_'+col] = multipletests(df[col], method="fdr_bh")[1]
        return(df)
    


    def list_most_significant_positions (self, n=10):
        pass

    def list_most_significant_references (self, n=10):
        pass

    ############################################################################

    def plot_signal (self, ref_id, start, end, figsize=(30,10), palette="Set2", plot_style="ggplot"):
        """
        Plot the dwell time and median intensity distribution position per positon in a split violin plot representation.
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate (Must be higher or equal to 0)
        end: End coordinate (included) (must be lower or equal to the reference length)
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style
            . See https://matplotlib.org/users/style_sheets.html
        """

        try:
            ref_fasta = self._fasta [ref_id]
        except KeyError:
            raise NanocomporeError ("Reference id not present in result database")
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
        if df.empty:
            raise NanocomporeError ("No data available for the selected interval")

        # Define ploting style
        with pl.style.context (plot_style):
            # Plot dwell and median
            fig, axes = pl.subplots(2, 1, figsize=figsize)
             ######### Add model value for median ?
            _ = sns.violinplot (x="pos", y="median", hue="sample", data=df, split=True, ax=axes[0], inner="quartile", bw=0.75, linewidth=1, palette=palette)
            _ = sns.violinplot (x="pos", y="dwell", hue="sample", data=df, split=True, ax=axes[1], inner="quartile", bw=0.75, linewidth=1, palette=palette)

            # Adjust display
            _ = axes[0].set_title (ref_id)
            _ = axes[0].get_xaxis().set_visible(False)
            _ = axes[0].set_ylabel ("Median Intensity")
            _ = axes[1].set_xticklabels (x_lab)
            _ = axes[1].set_ylabel ("Dwell Time")
            _ = axes[1].set_xlabel ("Reference position")

            pl.tight_layout()
            return (fig, axes)

    def plot_pvalue (self, ref_id, start=None, end=None, threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot"):
        """
        <Plot pvalues per position (by default plot all fields starting by "pvalue")
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        """
        try:
            ref_fasta = self._fasta [ref_id]
        except KeyError:
            raise NanocomporeError ("Reference id not present in result database")

        # Define start, end if not given
        if not start:
            start = 0
        if not end:
            end = len (ref_fasta)
        # Check start end
        if start > end:
            raise NanocomporeError ("End coordinate has to be higher or equal to start")
        if start < 0:
            raise NanocomporeError ("Coordinates have to be higher that 0")
        if end > len(ref_fasta):
            raise NanocomporeError ("Coordinates have to be lower than the ref_id sequence length ({})".format(len(ref_fasta)))

        # Parse line position per position
        d = OrderedDict ()
        ref_pos_dict = self [ref_id]
        for pos in range (start, end+1):
            # Collect results for position
            res_dict = OrderedDict ()
            #res_dict ["pos"] = pos
            if pos in ref_pos_dict:
                for k,v in ref_pos_dict[pos].items():
                    if k.startswith ("pvalue"): # Get every fields starting with "pvalue"
                        res_dict [k] = v
            d[pos] = res_dict

        # Cast collected results to dataframe
        df = pd.DataFrame.from_dict(d, orient="index")
        if df.empty:
            raise NanocomporeError ("No data available for the selected interval")

        # filling missing values and log transform the data
        df.fillna(1, inplace=True)
        df = -np.log(df)

        # Define ploting style
        with pl.style.context (plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            _ = sns.lineplot(data=df, palette=palette, ax=ax)
            _ = ax.axhline (y=-np.log(threshold), color="grey", linestyle=":", label="pvalue={}".format(threshold))
            _ = ax.legend ()
            _ = ax.set_ylabel ("-log (pvalue)")
            _ = ax.set_xlabel ("Reference position")
            _ = ax.set_title (ref_id)

            pl.tight_layout()
            return (fig, ax)
