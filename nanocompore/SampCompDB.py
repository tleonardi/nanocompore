# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, namedtuple
import shelve
from math import log
import re

# Third party
from pyfaidx import Fasta
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
import seaborn as sns
from bedparse import bedline
from statsmodels.stats.multitest import multipletests

# Local package
from nanocompore.common import counter_to_str, access_file, NanocomporeError
from nanocompore import models

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampCompDB (object):
    """ Wrapper over the result shelve SampComp """

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#
    def __init__(self, db_fn, fasta_fn, run_type="RNA"):
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

        if run_type == "RNA":
            self._model_dict = models.RNA_model_dict
        else:
            raise NanocomporeError ("Only RNA is implemented at the moment")

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
    def save_to_bed (self, output_fn, bedgraph=False, pvalue_field=None, pvalue_thr=0.01, sequence_context=0, convert=None, assembly=None, title=None):
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
            if title is not None:
                if not bedgraph:
                    bed_file.write('track type=bed name="%s" description="%s"\n'%(title,title))
                else:
                    bed_file.write('track type=bedGraph name="%s" description="%s"\n'%(title,title))
            for record in self.results[['chr', 'genomicPos', 'ref','strand']+[pvalue_field]].values.tolist():
                if not bedgraph and record[-1]<=pvalue_thr:
                    line=bedline([record[0], record[1], record[1]+sequence_context+1, record[2], -log(record[-1], 10), record[3]])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.name, line.score, line.strand))
                elif bedgraph:
                    line=bedline([record[0], record[1], record[1]+sequence_context+1, record[2], -log(record[-1], 10), record[3]])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.score))

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
            elif method in ["GMM"]:
                tests.append("pvalue_gmm")
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

    def plot_pvalue (self, ref_id, start=None, end=None, threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot", method=None):
        """
        Plot pvalues per position (by default plot all fields starting by "pvalue")
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        method: Limit the pvalue methods shown in the plot. Either a list of methods or a regular expression as a string.
        """
        try:
            ref_fasta = self._fasta [ref_id]
        except KeyError:
            raise NanocomporeError("Reference id not present in result database")

        # Define start, end if not given
        if not start:
            start = 0
        if not end:
            end = len (ref_fasta)
        # Check start end
        if start > end:
            raise NanocomporeError("End coordinate has to be higher or equal to start")
        if start < 0:
            raise NanocomporeError("Coordinates have to be higher that 0")
        if end > len(ref_fasta):
            raise NanocomporeError("Coordinates have to be lower than the ref_id sequence length ({})".format(len(ref_fasta)))

        try:
            ref_pos_dict = self.results.query('ref==@ref_id').set_index('pos').to_dict('index')
        except NameError:
            raise NanocomporeError("In order to plot adjusted pvalues you have to call the calculate_results() function first")

        # Make a list with all methods available
        methods=list(self.results)

        # If method not provided, set it to a default regex
        if method is None:
            method="adjusted_pvalue*"
        # Parse the method as regex if string
        if isinstance(method, str):
            r = re.compile(method)
            method = list(filter(r.match, methods))
        elif not isinstance(method, list):
            raise NanocomporeError("Method must be either a string or a list")
        
        for m in method:
            if m not in methods:
                raise NanocomporeError("Method %s is not the results dataframe"%m)

        # Parse line position per position
        d = OrderedDict ()

        for pos in range (start, end+1):
            # Collect results for position
            res_dict = OrderedDict ()
            if pos in ref_pos_dict:
                for k,v in ref_pos_dict[pos].items():
                    if k in method:
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

    def plot_signal (self, ref_id, start=None, end=None, split_samples=False, feature="intensity", figsize=(30,7), palette="Set2", plot_style="ggplot", bw=0.25):
        """
        Plot the dwell time and median intensity distribution position per positon in a split violin plot representation.
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate (Must be higher or equal to 0)
        end: End coordinate (included) (must be lower or equal to the reference length)
        split_samples: If samples for a same condition are represented separatly. If false they are merged per condition
        feature: Choose between "intensity" and "dwell"
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style
            . See https://matplotlib.org/users/style_sheets.html
        bw: Scale factor to use when computing the kernel bandwidth
        """

        # Get data
        ref_data, ref_fasta, start, end = self.__get_plot_data (ref_id, start, end)

        # Parse line position per position
        l = []
        valid=0
        x_ticks_list = []
        model_means_list = []

        if not feature in ["intensity", "dwell"]:
            raise NanocomporeError ("Feature either 'intensity' or 'dwell'")

        # Extract data from database if position in db
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                valid+=1
                for k1, v1 in ref_data[pos].items():
                    # Collect kmer seq
                    if k1 == "ref_kmer":
                        x_ticks_list.append("{}\n{}".format(pos, v1))
                        if feature == "intensity":
                            model_means_list.append (self._model_dict[v1][0])

                    # Collect dwell or median data
                    else:
                        for k2, v2 in v1.items():
                            for v3 in v2[feature]:
                                if split_samples:
                                    l.append ((pos, "{}_{}".format(k1, k2), v3))
                                else:
                                    l.append ((pos, k1, v3))
            else:
                l.append ((pos, None, None))
                x_ticks_list.append(str(pos))
                model_means_list.append(None)

        # Check that we found valid position and cast collected results to dataframe
        if not valid:
            raise NanocomporeError ("No data available for selected coordinates")
        df = pd.DataFrame (l, columns=["pos", "lab", feature])

        # Define ploting style
        with pl.style.context (plot_style):
            # Plot dwell and median
            fig, ax = pl.subplots (figsize=figsize)
            _ = sns.violinplot (
                x="pos",
                y=feature,
                hue="lab",
                data=df,
                split=not split_samples,
                ax=ax,
                inner="quartile",
                bw=bw,
                linewidth=1,
                palette=palette,
                scale="area")

            if feature == "intensity":
                _ = ax.plot (model_means_list, color="black", marker="x", label="Model Mean", linestyle="")
                _ = ax.set_ylabel ("Mean Intensity")
            elif feature == "dwell":
                _ = ax.set_ylabel ("Dwell Time")

            # Adjust display
            _ = ax.set_xlim (-1, end-start+1)
            _ = ax.set_xticklabels (x_ticks_list)
            _ = ax.set_title ("Reference:{}  Start:{}  End:{}".format(ref_id, start, end))
            _ = ax.set_xlabel ("Reference position")

            _ = ax.legend ()

            pl.tight_layout()
            return (fig, ax)

    def plot_coverage (self, ref_id, start=None, end=None, figsize=(30,5),  palette="Set2", plot_style="ggplot"):
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
        l=[]
        valid=0
        # Extract data from database if position in db
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                valid+=1
                for k1, v1 in ref_data[pos].items():
                    if k1 != "ref_kmer":
                        for k2, v2 in v1.items():
                            l.append ((pos, "{}_{}".format(k1, k2), v2["coverage"]))
            else:
                l.append ((pos, None, None))

        # Check that we found valid position and cast collected results to dataframe
        if not valid:
            raise NanocomporeError ("No data available for selected coordinates")
        df = pd.DataFrame (l, columns=["pos", "Sample", "cov"])

        # Define plotting style
        with pl.style.context (plot_style):
            fig, ax = pl.subplots (figsize=figsize)
            _ = sns.lineplot (
                x="pos",
                y="cov",
                hue="Sample",
                data=df,
                ax=ax,
                palette=palette,
                drawstyle="steps")

            _ = ax.set_ylim (0, None)
            _ = ax.set_xlim (start, end)
            _ = ax.set_title ("Reference:{}  Start:{}  End:{}".format(ref_id, start, end))
            _ = ax.set_ylabel ("Coverage")
            _ = ax.set_xlabel ("Reference position")
            pl.tight_layout()
            return (fig, ax)

    def plot_position(self, ref_id, pos=None, figsize=(30,10), colors=["dodgerblue", "salmon"], plot_style="ggplot"):
        """
        Plot the dwell time and median intensity at the given position as a scatter plot.
        ref_id: Valid reference id name in the database
        pos: Position of interest
        figsize: length and heigh of the output plot. Default=(30,10)
        colors: List of colors
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style
            . See https://matplotlib.org/users/style_sheets.html
        """

        # Get data
        ref_data, ref_fasta, start, end = self.__get_plot_data (ref_id, pos, pos)

        if pos not in ref_data.keys():
            raise NanocomporeError ("No data available for the selected position")

        # Parse line position per position
        lt = namedtuple ("lt", ["pos", "sample", "median", "dwell"])
        l = []
        for median, dwell in zip (ref_data[pos]["S1_median"], ref_data[pos]["S1_dwell"]):
            l.append (lt (pos, "S1", median, dwell))
        for median, dwell in zip (ref_data[pos]["S2_median"], ref_data[pos]["S2_dwell"]):
            l.append (lt (pos, "S2", median, dwell))

        # Check that we found valid position and cast collected results to dataframe
        df = pd.DataFrame (l)

        # Create x label including the original sequence and its position
        base=ref_fasta[pos]
        x_lab = "{}\n{}".format(pos, base)
        df['dwell'] = np.log10(df['dwell'])
        s1_median = df[df['sample']=="S1"]['median']
        s1_dwell = df[df['sample']=="S1"]['dwell']
        s2_median = df[df['sample']=="S2"]['median']
        s2_dwell = df[df['sample']=="S2"]['dwell']
        # Define ploting style
        with pl.style.context(plot_style):
            # Plot dwell and median
            fig, ax = pl.subplots(figsize=figsize)
            cmap1 = sns.light_palette(colors[0], as_cmap=True)
            cmap2 = sns.light_palette(colors[1], as_cmap=True)
            ax = sns.kdeplot(s1_median, s1_dwell, cmap=cmap1, label="S1", shade_lowest=False, legend=True)
            ax = sns.kdeplot(s2_median, s2_dwell, cmap=cmap2, label="S2", shade_lowest=False, legend=True)
            # Adjust display
            _ = ax.set_title ("%s\n%s%s"%(ref_id,base,pos))
            _ = ax.set_ylabel ("Dwell")
            _ = ax.set_xlabel ("Median Intensity")
            red_patch = mpatches.Patch(color=colors[0], label='S1')
            blue_patch = mpatches.Patch(color=colors[1], label='S2')
            pl.legend(handles=[red_patch, blue_patch])
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
