# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, namedtuple
import shelve
from math import log

# Third party
from pyfaidx import Fasta
import pandas as pd
import numpy as np
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

        # Try to get ref_id list and metadata from shelve db
        try:
            with shelve.open (db_fn, flag='r') as db:

                # Try to get metedata from db
                try:
                    metadata = db['__metadata']
                    self._comparison_method=metadata['comparison_method']
                    self._sequence_context=metadata['sequence_context']
                except KeyError:
                    raise NanocomporeError("The result database does not contain metadata")

                # Try to load read_ids
                self.ref_id_list = [k for k in db.keys() if k!='__metadata']
                if not self.ref_id_list:
                    raise NanocomporeError("The result database is empty")
        except:
            raise NanocomporeError("The result database cannot be opened")
        self._db_fn = db_fn

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
        return len (self.ref_id_list)-1

    def __iter__ (self):
        with shelve.open (self._db_fn, flag = "r") as db:
            for k, v in db.items():
                if not k == '__metadata':
                    yield (k, v)

    def __getitem__(self, items):
        with shelve.open (self._db_fn, flag = "r") as db:
            if items in db:
                return db[items]
            else:
                raise KeyError ("Item not found in the database")

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def save_to_bed (self, output_fn, bedgraph=False, pvalue_field=None, pvalue_thr=0.01, sequence_context=0, convert=None, assembly=None):
        """Saves the results object to BED6 format.
            bedgraph: save file in bedgraph format instead of bed
            pvalue_field: specifies what column to use as BED score (field 5, as -log10)
            pvalue_thr: only report positions with pvalue<=thr
            sequence_context: produce BED files for the given context
            convert: one of 'ensembl_to_ucsc' or 'ucsc_to_ensembl". Convert chromosome named between Ensembl and Ucsc conventions
            assembly: required if convert is used. One of "hg38" or "mm10"
        """
        if sequence_context != 0:
            pvalue_field=pvalue_field+"_context="+str(sequence_context)
        if pvalue_field not in self.results:
            raise NanocomporeError(("The field '%s' is not in the results" % pvalue_field))
        if "results" not in self.__dict__:
            raise NanocomporeError("Run calculate_results() before trying to save to bed")
        if convert not in [None, "ensembl_to_ucsc", "ucsc_to_ensembl"]:
            raise NanocomporeError("Convert value not valid")
        if convert is not None and assembly is None:
            raise NanocomporeError("The assembly argument is required in order to do the conversion. Choose one of 'hg38' or 'mm10' ")

        with open(output_fn, "w") as bed_file:
            for record in self.results[['chr', 'genomicPos', 'ref','strand']+[pvalue_field]].values.tolist():
                if not bedgraph and record[-1]<=pvalue_thr:
                    line=bedline([record[0], record[1], record[1]+sequence_context+1, record[2], -log(record[-1], 10), record[3]])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s %s %s %s %s %s\n" % (line.chr, line.start, line.end, line.name, line.score, line.strand))
                elif bedgraph:
                    line=bedline([record[0], record[1], record[1]+sequence_context+1, record[2], -log(record[-1], 10), record[3]])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s %s %s %s\n" % (line.chr, line.start, line.end, line.score))

    def calculate_results(self, bed_fn=None, adjust=True, methods=None):
        # Compose a lists with the name of the results
        tests=[]
        if methods is None:
            methods = self._comparison_method
        for method in methods:
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

        # We open the DB rather that calling __getitem__ to avoid
        # opening and closing the file for every transcript
        with shelve.open(self._db_fn, flag = "r") as db:
            df = pd.DataFrame([dict({a:b for a,b in v.items() if a in tests}, ref=ref_id, pos=k)  for ref_id, rec in db.items() for k,v in rec.items() if ref_id!="__metadata" ]).fillna(1)

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
            # This is very inefficient. We should get chr and strand only once per transcript, ideally when writing the BED file
            df['chr'] = df.apply(lambda row: bed_annot[row['ref']].chr,axis=1)
            df['strand'] = df.apply(lambda row: bed_annot[row['ref']].strand,axis=1)
            df=df[['ref', 'pos', 'chr', 'strand', 'genomicPos']+tests]

        if adjust:
            for col in tests:
                df['adjusted_'+col] = multipletests(df[col], method="fdr_bh")[1]
        self.results = df

    def list_most_significant_positions (self, n=10):
        pass

    def list_most_significant_references (self, n=10):
        pass

    #~~~~~~~~~~~~~~PLOTING METHODS~~~~~~~~~~~~~~#

    def plot_pvalue (self, ref_id, start=None, end=None, adjusted_pvalues=True, threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot"):
        """
        Plot pvalues per position (by default plot all fields starting by "pvalue")
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        adjusted_pvalues: plot adjusted pvalues. Requires a results slot to be defined.
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
        if adjusted_pvalues:
            try:
                ref_pos_dict = self.results.query('ref==@ref_id').set_index('pos').to_dict('index')
                pvalue_selector="adjusted_"
            except NameError:
                raise NanocomporeError("In order to plot adjusted pvalues you have to call the results() function first")
        else:
                ref_pos_dict = self[ref_id]
                pvalue_selector="pvalue_"

        for pos in range (start, end+1):
            # Collect results for position
            res_dict = OrderedDict ()
            if pos in ref_pos_dict:
                for k,v in ref_pos_dict[pos].items():
                    if k.startswith (pvalue_selector): # Get every fields starting with "pvalue"
                        res_dict [k] = v
            d[pos] = res_dict

        # Create x label including the original sequence and its position
        x_lab = []
        for pos, base in zip (range (start, end+1), ref_fasta[start:end+1]):
            x_lab.append ("{}\n{}".format(pos, base))

        # Cast collected results to dataframe
        df = pd.DataFrame.from_dict(d, orient="index")
        if df.empty:
            raise NanocomporeError ("No data available for the selected interval")

        # filling missing values and log transform the data
        df.fillna(1, inplace=True)
        df = -np.log10(df)

        # Define plotting style
        with pl.style.context (plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            _ = sns.lineplot(data=df, palette=palette, ax=ax, dashes=False)
            _ = ax.axhline (y=-np.log10(threshold), color="grey", linestyle=":", label="pvalue={}".format(threshold))
            _ = ax.legend ()
            _ = ax.set_ylabel ("-log (pvalue)")
            _ = ax.set_xlabel ("Reference position")
            _ = ax.set_title (ref_id)
            _ = ax.set_xlim (start, end)
            if end-start<30:
                _ = ax.set_xticks(df.index)
                _ = ax.set_xticklabels(x_lab)
            pl.tight_layout()
            return (fig, ax)

    def plot_signal (self, ref_id, start=None, end=None, figsize=(30,10), colors=["dodgerblue", "salmon"], plot_style="ggplot"):
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

        # Get data
        ref_data, ref_fasta, start, end = self.__get_plot_data (ref_id, start, end)

        # Parse line position per position
        lt = namedtuple ("lt", ["pos", "sample", "median", "dwell"])
        l = []
        valid_pos=0
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                valid_pos+=1
                for median, dwell in zip (ref_data[pos]["S1_median"], ref_data[pos]["S1_dwell"]):
                    l.append (lt (pos, "S1", median, dwell))
                for median, dwell in zip (ref_data[pos]["S2_median"], ref_data[pos]["S2_dwell"]):
                    l.append (lt (pos, "S2", median, dwell))

            # If not coverage for position, just fill in with empty values
            else:
                for sample in ("S1", "S2"):
                    l.append (lt (pos, sample, None, None))

        # Check that we found valid position and cast collected results to dataframe
        if not valid_pos:
            raise NanocomporeError ("No data available for the selected interval")
        df = pd.DataFrame (l)

        # Create x label including the original sequence and its position
        x_lab = []
        for pos, base in zip (range (start, end+1), ref_fasta[start:end+1]):
            x_lab.append ("{}\n{}".format(pos, base))

        # Define ploting style
        with pl.style.context (plot_style):
            # Plot dwell and median
            fig, axes = pl.subplots(2, 1, figsize=figsize)
            _ = sns.violinplot (x="pos", y="median", hue="sample", data=df, split=True, ax=axes[0], inner="quartile", bw=0.75, linewidth=1, palette=colors)
            _ = sns.violinplot (x="pos", y="dwell", hue="sample", data=df, split=True, ax=axes[1], inner="quartile", bw=0.75, linewidth=1, palette=colors)

            # Adjust display
            _ = axes[0].set_title (ref_id)
            _ = axes[0].get_xaxis().set_visible(False)
            _ = axes[0].set_ylabel ("Median Intensity")
            _ = axes[1].set_xticklabels (x_lab)
            _ = axes[1].set_ylabel ("Dwell Time")
            _ = axes[1].set_xlabel ("Reference position")

            pl.tight_layout()
            return (fig, axes)

    def plot_coverage (self, ref_id, start=None, end=None, figsize=(30,5), colors=["dodgerblue", "salmon"], plot_style="ggplot"):
        """
        Plot pvalues per position (by default plot all fields starting by "pvalue")
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        figsize: length and heigh of the output plot. Default=(30,10)
        colors: list of 2 colors
            see https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        """
        # Get data
        ref_data, ref_fasta, start, end = self.__get_plot_data (ref_id, start, end)


        # Parse line position per position
        lt = namedtuple ("lt", ["pos", "S1", "S2"])
        l = []
        valid_pos = 0
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                l.append (lt (pos, ref_data[pos]["S1_coverage"], ref_data[pos]["S2_coverage"]))
                valid_pos+=1
            else:
                l.append (lt (pos, None, None))

        # Check that we found valid position and cast collected results to dataframe
        if not valid_pos:
            raise NanocomporeError ("No data available for the selected interval")
        df = pd.DataFrame (l)

        # Define plotting style
        with pl.style.context (plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            _ = ax.plot (df["pos"], df["S1"], label="S1", color=colors[0])
            _ = ax.plot (df["pos"], df["S2"], label="S2", color=colors[1])
            _ = ax.set_ylabel ("Coverage")
            _ = ax.set_xlabel ("Reference position")
            _ = ax.set_xlim (start, end)
            _ = ax.set_title (ref_id)
            _ = ax.legend ()
            pl.tight_layout()
            return (fig, ax)

    #~~~~~~~~~~~~~~PRIVATE  METHODS~~~~~~~~~~~~~~#
    def __get_plot_data (self, ref_id, start, end):
        """
        Private function to verify and extract info for the ploting functions
        """

        # Extract data for ref_id
        try:
            ref_data = self[ref_id]
        except KeyError:
            raise NanocomporeError ("Reference id not present in result database")
        # Get corresponding fasta record
        try:
            ref_fasta = self._fasta [ref_id]
        except KeyError:
            raise NanocomporeError ("Reference id not present in fasta reference database")
        # Get start and end
        if not start:
            start = 0
        if not end:
            end = len (ref_fasta)
        if start > end:
            raise NanocomporeError ("End coordinate has to be higher or equal to start")
        if start < 0:
            raise NanocomporeError ("Coordinates have to be higher that 0")
        if end > len(ref_fasta):
            raise NanocomporeError ("Coordinates have to be lower than the ref_id sequence length ({})".format(len(ref_fasta)))

        return (ref_data, ref_fasta, start, end)
