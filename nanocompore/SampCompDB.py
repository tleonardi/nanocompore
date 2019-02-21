# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import shelve
from dbm import error as dbm_error
from math import log
import re
import sys
from pkg_resources import resource_filename

# Third party
from pyfaidx import Fasta
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
import seaborn as sns
from bedparse import bedline
from statsmodels.stats.multitest import multipletests
from sklearn.mixture.gaussian_mixture import GaussianMixture
from sklearn.preprocessing import scale as scale

# Local package
from nanocompore.common import *

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampCompDB(object):
    """ Wrapper over the result shelve SampComp """

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#
    def __init__(self, db_fn, fasta_fn, bed_fn=None, run_type="RNA"):
        """
        Import a shelve db and a fasta reference file. Automatically returned by SampComp
        Can also be manually created from an existing shelve db output
        db_fn: Path where to write the result database
        fasta_fn: Path to a fasta file corresponding to the reference used for read alignemnt
        bed_fn: Path to a BED file containing the annotation of the transcriptome used as reference when mapping
        """

        # Try to get ref_id list and metadata from shelve db
        try:
            with shelve.open(db_fn, flag='r') as db:
                # Try to get metadata from db
                try:
                    metadata = db['__metadata']
                    self._comparison_method = metadata['comparison_method']
                    self._sequence_context = metadata['sequence_context']
                    self._min_coverage = metadata['min_coverage']
                    self._n_samples = metadata['n_samples']
                    self._pvalue_tests = metadata['pvalue_tests']
                except KeyError:
                    raise NanocomporeError("The result database does not contain metadata")
                # Try to load read_ids
                self.ref_id_list = [k for k in db.keys() if k!='__metadata']
                if not self.ref_id_list:
                    raise NanocomporeError("The result database is empty")
                # finally save db path
                self._db_fn = db_fn

        except dbm_error:
            raise NanocomporeError("The result database cannot be opened")

        # Test is Fasta can be opened
        try:
            with Fasta(fasta_fn):
                self._fasta_fn = fasta_fn
        except IOError:
            raise NanocomporeError("The fasta file cannot be opened")

        # Define model depending on run_type
        if run_type == "RNA":
            model_fn = resource_filename("nanocompore", "models/kmers_model_RNA_r9.4_180mv.tsv")
            self._model_df = pd.read_csv(model_fn, sep="\t", comment="#", index_col=0)
        else:
            raise NanocomporeError("Only RNA is implemented at the moment")

        self.bed_fn = bed_fn

        #Create results df with adjusted p-values
        if self._pvalue_tests:
            self.results = self.__calculate_results(adjust=True)

    def __repr__(self):
        """readable description of the object"""
        return "[{}] Number of references: {}\n".format(self.__class__.__name__, len(self))

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __len__(self):
        return len(self.ref_id_list)-1

    def __iter__(self):
        with shelve.open(self._db_fn, flag = "r") as db:
            for k, v in db.items():
                if not k == '__metadata':
                    yield(k, v)

    def __getitem__(self, items):
        with shelve.open(self._db_fn, flag = "r") as db:
            if items in db:
                return db[items]
            else:
                raise KeyError("Item not found in the database")

    #~~~~~~~~~~~~~~PRIVATE  METHODS~~~~~~~~~~~~~~#

    def __calculate_results(self, adjust=True):
        """"""

        # Collect all pvalue results in a dataframe
        l = []
        for ref_id, ref_pos_list in self:
            for pos, pos_dict in enumerate(ref_pos_list):
                if "txComp" in pos_dict:
                    row_dict = OrderedDict()
                    row_dict["ref_id"] = ref_id
                    row_dict["pos"] = pos
                    row_dict["ref_kmer"] = pos_dict["ref_kmer"]
                    for test in self._pvalue_tests:
                        if test in pos_dict["txComp"]:
                            pval = pos_dict["txComp"][test]
                            row_dict[test] = pval
                    l.append(row_dict)
        df = pd.DataFrame(l)

        if self.bed_fn:
            bed_annot={}
            try:
                with open(self.bed_fn) as tsvfile:
                    for line in tsvfile:
                        record_name=line.split('\t')[3]
                        if( record_name in self.ref_id_list):
                            bed_annot[record_name]=bedline(line.split('\t'))
            except:
                raise NanocomporeError("Can't open BED file")
            if len(bed_annot) != len(self.ref_id_list):
                raise NanocomporeError("Some references are missing from the BED file provided")

            df['genomicPos'] = df.apply(lambda row: bed_annot[row['ref_id']].tx2genome(coord=row['pos'], stranded=True),axis=1)
            # This is very inefficient. We should get chr and strand only once per transcript, ideally when writing the BED file
            df['chr'] = df.apply(lambda row: bed_annot[row['ref_id']].chr,axis=1)
            df['strand'] = df.apply(lambda row: bed_annot[row['ref_id']].strand,axis=1)
            df=df[['ref_id', 'pos', 'chr', 'strand', 'genomicPos', 'ref_kmer']+self._pvalue_tests]
        else:
            df=df[['ref_id', 'pos', 'ref_kmer']+self._pvalue_tests]

        if adjust:
            for col in self._pvalue_tests:
                df[col] = self.__multipletests_filter_nan(df[col], method="fdr_bh")
        return df

    def __get_kmer_list(self, ref_id, start, end, kmer_size=5):
        """ Extract fasta record corresponding to ref with error handling """
        try:
            with Fasta(self._fasta_fn) as fasta:
                fasta =(fasta [ref_id])
                seq = str(fasta[start:end+5])
                kmer_list = []
                for i in range(end-start):
                    kmer_list.append(seq[i:i+5])
                return kmer_list
        except KeyError:
            raise NanocomporeError("Reference id not present in fasta file")

    def __get_positions(self, ref_id, start=None, end=None):
        """ Verify start and end and if not available try to infer them"""
        try:
            with Fasta(self._fasta_fn) as fasta:
                max_len = len(fasta [ref_id])-4
        except KeyError:
            raise NanocomporeError("Reference id not present in fasta file")
        if not start or start < 0:
            start = 0
        if not end or end > max_len:
            end = max_len
        if start > end:
            raise NanocomporeError("End coordinate has to be higher or equal to start")
        return(start, end)

    @staticmethod
    def __color_generator(palette, n):
        pal = sns.mpl_palette(palette, n)
        for i in range(n):
            yield(pal[i])

    @staticmethod
    def __multipletests_filter_nan(pvalues, method="fdr_bh"):
        """
        Performs p-value correction for multiple hypothesis testing
        using the method specified. The pvalues list can contain
        np.nan values, which are ignored during p-value correction.
        test: input=[0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan]
        out: array([0.13333333, 0.016     ,        nan, 0.016     , 0.5       ,
        0.45714286, 0.016     , 0.008     ,        nan,        nan,
        0.016     ,        nan])
        """
        pvalues_no_nan = [p for p in pvalues if not np.isnan(p)]
        corrected_p_values = multipletests(pvalues_no_nan, method=method)[1]
        for i, p in enumerate(pvalues):
            if np.isnan(p):
                corrected_p_values=np.insert(corrected_p_values, i, np.nan, axis=0)
        return(corrected_p_values)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def save_to_bed(self, output_fn, bedgraph=False, pvalue_field=None, pvalue_thr=0.01, span=5, convert=None, assembly=None, title=None):
        """Saves the results object to BED6 format.
            bedgraph: save file in bedgraph format instead of bed
            pvalue_field: specifies what column to use as BED score (field 5, as -log10)
            pvalue_thr: only report positions with pvalue<=thr
            span: The size of each BED feature.
                  If size=5 (default) features correspond to kmers. If size=1 features correspond to the first base of each kmer.
            convert: one of 'ensembl_to_ucsc' or 'ucsc_to_ensembl". Convert chromosome named between Ensembl and Ucsc conventions
            assembly: required if convert is used. One of "hg38" or "mm10"
        """
        if self.bed_fn is None:
            raise NanocomporeError("In order to generate a BED file SampCompDB needs to be initialised with a transcriptome BED")
        if span < 1:
            raise NanocomporeError("span has to be >=1")
        if span != 5 and bedgraph:
            raise NanocomporeError("Span is ignored when generating bedGraph files")
        if pvalue_field not in self.results:
            raise NanocomporeError(("The field '%s' is not in the results" % pvalue_field))
        if "results" not in self.__dict__:
            raise NanocomporeError("It looks like there's not results slot in SampCompDB")
        if convert not in [None, "ensembl_to_ucsc", "ucsc_to_ensembl"]:
            raise NanocomporeError("Convert value not valid")
        if convert is not None and assembly is None:
            raise NanocomporeError("The assembly argument is required in order to do the conversion. Choose one of 'hg38' or 'mm10' ")

        with open(output_fn, "w") as bed_file:
            if title is not None:
                if not bedgraph:
                    bed_file.write('track type=bed name="%s" description="%s"\n'%(title,pvalue_field))
                else:
                    bed_file.write('track type=bedGraph name="%s" description="%s"\n'%(title,pvalue_field))

            Record = namedtuple('Record', ['chr', 'genomicPos', 'ref_id','strand', 'ref_kmer', pvalue_field ])
            for record in self.results[ list(Record._fields) ].itertuples(index=False, name="Record"):
                pvalue = getattr(record, pvalue_field)
                if np.isnan(pvalue):
                    pvalue=0
                elif pvalue < sys.float_info.min:
                    pvalue = -log(sys.float_info.min, 10)
                else:
                    pvalue=-log(pvalue, 10)
                if not bedgraph and pvalue >= -log(pvalue_thr, 10):
                    if record.strand == "+":
                        line=bedline([record.chr, record.genomicPos, record.genomicPos+span, f"{record.ref_id}_{record.ref_kmer}", pvalue, record.strand])
                    else:
                        line=bedline([record.chr, record.genomicPos-(span-1), record.genomicPos+1, f"{record.ref_id}_{record.ref_kmer}", pvalue, record.strand])

                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.name, line.score, line.strand))
                elif bedgraph:
                    if record.strand == "+":
                        line=bedline([record.chr, record.genomicPos+2, record.genomicPos+3, f"{record.ref_id}_{record.ref_kmer}", pvalue, record.strand])
                    else:
                        line=bedline([record.chr, record.genomicPos-2, record.genomicPos-1, f"{record.ref_id}_{record.ref_kmer}", pvalue, record.strand])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.score))

    def save_report(self, output_fn=None):
        """Saves an extended tabular report
        """
        if output_fn is None:
            fp = sys.stdout
        elif isinstance(output_fn, str):
            try:
                fp = open(output_fn, "w")
            except:
                raise NanocomporeError("Error opening output file %s"%output_fn)
        else:
            raise NanocomporeError("output_fn needs to be a string or None")

        headers = ['chr', 'pos', 'ref_id','strand', 'ref_kmer']+self._pvalue_tests
        # Read extra GMM info from the shelve
        if "GMM" in self._comparison_method:
            headers += ['LOR', 'clusters']
            gmm_info=OrderedDict()
            for tx, refpos in self:
                gmm_info[tx] = {k:{'lor': v['txComp']['GMM_model'][1], 'clusters':v['txComp']['GMM_model'][3]} for k,v in enumerate(refpos) if "txComp" in v and "GMM_model" in v['txComp']}
        fp.write('\t'.join([ str(i) for i in headers ])+'\n')
        for record in self.results[['chr', 'pos', 'ref_id','strand', 'ref_kmer']+self._pvalue_tests].values.tolist():
            if "GMM" in self._comparison_method:
                try:
                    lor = gmm_info[record[2]][record[1]]['lor']
                    clusters = gmm_info[record[2]][record[1]]['clusters']
                    record += [lor, clusters]
                except KeyError:
                    record += ["nan", "nan"]
            fp.write('\t'.join([ str(i) for i in record ])+'\n')
        fp.close()


    def list_most_significant_positions(self, n=10):
        pass

    def list_most_significant_references(self, n=10):
        pass

    #~~~~~~~~~~~~~~PLOTTING METHODS~~~~~~~~~~~~~~#
    def plot_pvalue( self, ref_id, start=None, end=None, kind="lineplot", threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot", tests=None):
        """
        Plot pvalues per position (by default plot all fields starting by "pvalue")
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        kind: kind of plot to represent the data (lineplot or barplot)
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        tests: Limit the pvalue methods shown in the plot. Either a list of methods or a string coresponding to a part of the name
        """
        # Extract fasta and positions
        start, end = self.__get_positions(ref_id, start, end)
        kmer_list = self.__get_kmer_list(ref_id, start, end)
        pos_list = list(range(start, end))

        # Extract data from result dataframe
        df = self.results.query("ref_id==@ref_id and @start<=pos<@end")
        if df.empty:
            raise NanocomporeError("No data available for the selected interval")

        # list tests methods
        if not tests:
            tests = self._pvalue_tests
        else:
            if isinstance(tests, str):
                tests = tests.split(",")
            elif not isinstance(tests, list):
                raise NanocomporeError("Method must be either a string or a list")
            tests_set = set()
            for test in tests:
                for available_test in self._pvalue_tests:
                    if test in available_test:
                        tests_set.add(available_test)
            if not tests_set:
                raise NanocomporeError("No matching tests found in the list of available tests: {}".format(self._pvalue_tests))
            tests = sorted(list(tests_set))

        # Collect data for interval in a numpy array
        array = np.zeros(end-start, dtype=[(test, np.float) for test in tests])
        for id, row in df.iterrows():
            for test in tests:
                array[row["pos"]-start][test] = -np.log10(row[test])
        # Cast collected results to dataframe
        df = pd.DataFrame(array, index=pos_list)

        # Define plotting style
        with pl.style.context(plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            if kind == "lineplot":
                _ = sns.lineplot(data=df, palette=palette, ax=ax, dashes=False)
                if end-start<30:
                    _ = ax.set_xticks(pos_list)
                    _ = ax.set_xticklabels([ "{}\n{}".format(i, j) for i,j in zip(pos_list, kmer_list)])

            elif kind == "barplot":
                df = df.reset_index()
                df = df.melt(id_vars="index", var_name="method", value_name="pvalue")
                _ = sns.barplot(x="index", y="pvalue", hue="method", data=df, palette=palette, ax=ax)
                if end-start<30:
                    _ = ax.set_xticklabels([ "{}\n{}".format(i, j) for i,j in zip(pos_list, kmer_list)])
                else:
                    breaks = list(range(start, end+1, (end-start)//10))
                    _ = ax.set_xticklabels([ i if i in breaks else "" for i in pos_list])

            _ = ax.axhline(y=-np.log10(threshold), color="grey", linestyle=":", label="pvalue={}".format(threshold))
            _ = ax.legend(bbox_to_anchor=(1, 1), loc=2, facecolor="white", frameon=False)
            _ = ax.set_ylabel("-log (pvalue)")
            _ = ax.set_xlabel("Reference position")
            _ = ax.set_title(ref_id)
            pl.tight_layout()
            return(fig, ax)

    def plot_signal(self, ref_id, start=None, end=None, kind="violinplot", split_samples=False, figsize=(30,10), palette="Set2", plot_style="ggplot"):
        """
        Plot the dwell time and median intensity distribution position per positon in a split violin plot representation.
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate (Must be higher or equal to 0)
        end: End coordinate (included) (must be lower or equal to the reference length)
        kind: Kind of plot. Can be violinplot, boxenplot or swarmplot
        split_samples: If samples for a same condition are represented separatly. If false they are merged per condition
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style
            . See https://matplotlib.org/users/style_sheets.html
        bw: Scale factor to use when computing the kernel bandwidth
        """

        # Extract data for ref_id
        ref_data = self[ref_id]
        start, end = self.__get_positions(ref_id, start, end)

        # Parse line position per position
        l_intensity = []
        l_dwell = []
        x_ticks_list = []
        model_intensity_list = []
        model_dwell_list = []

        # Extract data from database if position in db
        for pos in np.arange(start, end):
            ref_kmer=ref_data[pos]['ref_kmer']
            x_ticks_list.append("{}\n{}".format(pos, ref_kmer))
            model_intensity_list.append(self._model_df.loc[ref_kmer]["model_intensity_median"])
            model_dwell_list.append(np.log10(self._model_df.loc[ref_kmer]["model_dwell_median"]))

            for cond_lab, cond_dict in ref_data[pos]['data'].items():
                for sample_lab, sample_val in cond_dict.items():
                    lab = "{}_{}".format(cond_lab, sample_lab) if split_samples else cond_lab

                    # Add intensity and dwell values to list for curent pos / lab
                    if not sample_val["intensity"]:
                        l_intensity.append((pos, lab, None))
                    for value in sample_val["intensity"]:
                        l_intensity.append((pos, lab, value))
                    if not sample_val["dwell"]:
                        l_dwell.append((pos, lab, None))
                    for value in sample_val["dwell"]:
                        l_dwell.append((pos, lab, np.log10(value)))

        # Define ploting style
        with pl.style.context(plot_style):
            fig, (ax1, ax2) = pl.subplots(2,1, figsize=figsize, sharex=True)

            for ax, l, model in ((ax1,l_intensity,model_intensity_list), (ax2,l_dwell,model_dwell_list)):
                # Plot values
                df = pd.DataFrame(l, columns=["pos", "lab", "value"])
                if kind == "violinplot":
                    _ = sns.violinplot(x="pos", y="value", hue="lab", data=df, ax=ax, split=not split_samples, inner="quartile",
                                        bw=0.25, linewidth=1, scale="area", palette=palette, zorder=0)
                elif kind == "boxenplot":
                    _ = sns.boxenplot(x="pos", y="value", hue="lab", data=df, ax=ax, scale="area", palette=palette, zorder=0)
                elif kind == "swarmplot":
                    _ = sns.swarmplot(x="pos", y="value", hue="lab", data=df, ax=ax, dodge=True, palette=palette, zorder=0)
                else:
                    raise NanocomporeError("Not a valid plot kind {}".format(kind))
                # Plot model
                _ = ax.plot(model, color="black", marker="x", label="Model Median", linestyle="", zorder=1)

            # Adjust display
            _ = ax1.set_xlabel("")
            _ = ax2.set_xlabel("Reference position")
            _ = ax1.set_ylabel("Mean Intensity")
            _ = ax2.set_ylabel("log10 (Dwell Time)")
            _ = ax1.legend(bbox_to_anchor=(1, 1), loc=2, facecolor="white", frameon=False)
            _ = ax2.legend("")
            _ = ax2.set_xlim(-1, end-start)
            _ = ax2.set_xticklabels(x_ticks_list)
            _ = fig.suptitle("Reference:{}  Start:{}  End:{}".format(ref_id, start, end), y=1.01, fontsize=18)

            pl.tight_layout()
            return(fig, (ax1, ax2))

    def plot_coverage(self, ref_id, start=None, end=None, scale=False, split_samples=False, figsize=(30,5), palette="Set2", plot_style="ggplot"):
        """
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        """
        # Extract data for ref_id
        ref_data = self[ref_id]
        start, end = self.__get_positions(ref_id, start, end)

        # Parse data from database
        l = []
        for pos in np.arange(start, end):
            for cond_lab, cond_dict in ref_data[pos]['data'].items():
                for sample_lab, sample_val in cond_dict.items():
                    l.append((pos, "{}_{}".format(cond_lab, sample_lab), sample_val["coverage"]))

        # Cast collected results to dataframe
        df = pd.DataFrame(l, columns=["pos", "sample", "cov"])
        if scale:
            df['cov'] = df.groupby('sample')['cov'].apply(lambda x: x/max(x))

        # Define plotting style
        with pl.style.context(plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            _ = sns.lineplot( x="pos", y="cov", hue="sample", data=df, ax=ax, palette=palette, drawstyle="steps")
            if not scale:
                _ = ax.axhline(y=self._min_coverage, linestyle=":", color="grey", label="minimal coverage")

            _ = ax.set_ylim(0, None)
            _ = ax.set_xlim(start, end-1)
            _ = ax.set_title("Reference:{}  Start:{}  End:{}".format(ref_id, start, end), fontsize=18)
            _ = ax.set_ylabel("Coverage")
            _ = ax.set_xlabel("Reference position")
            _ = ax.legend(bbox_to_anchor=(1, 1), loc=2, facecolor="white", frameon=False)

            pl.tight_layout()
            return(fig, ax)

    def plot_bleeding_hulk(self, ref_id, start=None, end=None, split_samples=False, figsize=(30,10)):
        self.plot_kmers_stats(ref_id, start, end, split_samples, figsize, "Accent")

    def plot_kmers_stats(self, ref_id, start=None, end=None, split_samples=False, figsize=(30,10), palette="Accent", plot_style="ggplot"):
        """
        ref_id: Valid reference id name in the database
        start: Start coordinate. Default=0
        end: End coordinate (included). Default=reference length
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style. Default="ggplot"
            . See https://matplotlib.org/users/style_sheets.html
        """
        # Extract data for ref_id
        ref_data = self[ref_id]
        start, end = self.__get_positions(ref_id, start, end)

        # Parse data from database
        d = OrderedDict()
        for pos in np.arange(start, end):
            for cond_lab, cond_dict in ref_data[pos]["data"].items():
                for samp_lab, sample_val in cond_dict.items():
                    lab = "{}_{}".format(cond_lab, samp_lab) if split_samples else cond_lab

                    # Create dict arborescence
                    if not lab in d:
                        d[lab] = OrderedDict()
                    if not pos in d[lab]:
                        d[lab][pos] = {"valid":0,"NNNNN":0,"mismatching":0,"missing":0}

                    # Fill-in with values
                    d[lab][pos]["valid"] += sample_val["kmers_stats"]["valid"]
                    d[lab][pos]["NNNNN"] += sample_val["kmers_stats"]["NNNNN"]
                    d[lab][pos]["mismatching"] += sample_val["kmers_stats"]["mismatching"]
                    d[lab][pos]["missing"] += sample_val["kmers_stats"]["missing"]

        with pl.style.context(plot_style):
            fig, axes = pl.subplots(len(d),1, figsize=figsize)
            for ax, (lab, pos_dict) in zip(axes, d.items()):
                sample_df = pd.DataFrame.from_dict(pos_dict, orient="index")

                _ = sample_df.plot.area(ax=ax, colormap=palette, legend=False)
                _ = ax.set_title(lab)
                _ = ax.set_ylabel("Coverage")
                _ = ax.set_xlim(start, end-1)

            _ = axes[-1].set_xlabel("Reference position")
            _ = axes[0].legend(bbox_to_anchor=(1, 1), loc=2, facecolor="white", frameon=False)
            _ = fig.suptitle("Reference:{}  Start:{}  End:{}".format(ref_id, start, end), y=1.02, fontsize=18)
            pl.tight_layout()

        return(fig, axes)

    def plot_position(self, ref_id, pos=None, split_samples=False, figsize=(30,10), palette="Set2",  plot_style="ggplot", xlim=(None,None), ylim=(None,None), alpha=0.3, pointSize=20, scatter=True, kde=True, model=False, gmm_levels=50):
        """
        Plot the dwell time and median intensity at the given position as a scatter plot.
        ref_id: Valid reference id name in the database
        pos: Position of interest
        split_samples: If True, samples for a same condition are represented separately. If False, they are merged per condition
        figsize: length and heigh of the output plot. Default=(30,10)
        palette: Colormap. Default="Set2"
            see https://matplotlib.org/users/colormaps.html, https://matplotlib.org/examples/color/named_colors.html
        plot_style: Matplotlib plotting style.
            See https://matplotlib.org/users/style_sheets.html
        xlim: A tuple of explicit limits for the x axis
        ylim: A tuple of explicit limits for the y axis
        kde: plot the KDE of the intensity/dwell bivarariate distributions in the two samples
        scatter: if True, plot the individual data points
        pointSize: int specifying the point size for the scatter plot
        model: If true, plot the GMM density estimate
        gmm_levels: number of contour lines to use for the GMM countour plot
        """
        # Extract data for ref_id
        ref_data = self[ref_id]

        # Check that position is valid
        if not isinstance(pos, int):
            raise NanocomporeError("pos must be a single position")
        if pos > len(ref_data):
            raise NanocomporeError("Position out of range")
        # if not ref_data[pos]['data']["intensity"] or not ref_data[pos]['data']["dwell"]:
        #     raise NanocomporeError("No data found for selected position")

        # Extract data from database if position in db
        ref_kmer = ref_data[pos]['ref_kmer']
        data = ref_data[pos]['data']

        # Sample colors in palette
        col_gen = self.__color_generator(palette=palette, n=self._n_samples if split_samples else 2)

        # Collect and transform data in dict
        plot_data_dict = OrderedDict()
        for cond_lab, cond_dict in ref_data[pos]["data"].items():
            if split_samples:
                for samp_lab, sample_val in cond_dict.items():
                    plot_data_dict["{}_{}".format(cond_lab, samp_lab)] = {
                        "intensity":scale(sample_val["intensity"]),
                        "dwell":scale(np.log10(sample_val["dwell"])),
                        "color":next(col_gen)}
            else:
                intensity_list = []
                dwell_list = []
                for samp_lab, sample_val in cond_dict.items():
                    intensity_list.append(sample_val["intensity"])
                    dwell_list.append(sample_val["dwell"])
                plot_data_dict[cond_lab] = {
                    "intensity":scale(np.concatenate(intensity_list)),
                    "dwell":scale(np.log10(np.concatenate(dwell_list))),
                    "color":next(col_gen)}

        # Add GMM model if required and available
        if model and 'txComp' in ref_data[pos] and 'GMM_model' in ref_data[pos]['txComp'] and isinstance(model, GaussianMixture):
            model = ref_data[pos]['txComp']['GMM_model'][4]
            condition_labels = tuple(data.keys())
            global_intensity = scale(np.concatenate(([v['intensity'] for v in data[condition_labels[0]].values()]+[v['intensity'] for v in data[condition_labels[1]].values()]), axis=None))
            global_dwell = scale(np.log10(np.concatenate(([v['dwell'] for v in data[condition_labels[0]].values()]+[v['dwell'] for v in data[condition_labels[1]].values()]), axis=None)))
            x = np.linspace(min(global_intensity), max(global_intensity), num=1000)
            y = np.linspace(min(global_dwell), max(global_dwell), num=1000)
            X, Y = np.meshgrid(x, y)
            XX = np.array([X.ravel(), Y.ravel()]).T
            Z = -model.score_samples(XX)
            Z = Z.reshape(X.shape)
        else:
            model = None

        # plot collected data
        with pl.style.context(plot_style):
            fig, ax = pl.subplots(figsize=figsize)

            for label, d in plot_data_dict.items():
                if kde:
                    _ = sns.kdeplot(
                        data=d["intensity"],
                        data2=d["dwell"],
                        cmap=sns.light_palette(d["color"], as_cmap=True),
                        ax=ax,
                        clip=((min(d["intensity"]), max(d["intensity"])), (min(d["dwell"]),max(d["dwell"]))))
                if scatter:
                    _ = ax.scatter(
                        x=d["intensity"],
                        y=d["dwell"],
                        color=d["color"],
                        label=label,
                        alpha=alpha,
                        s=pointSize)
            if model:
                _ = ax.contour(X, Y, Z, levels=gmm_levels, alpha=alpha, colors="black")

            # Adjust display
            _ = ax.set_title("%s\n%s (%s)"%(ref_id,pos, ref_kmer))
            _ = ax.set_ylabel("log10 (Dwell Time)")
            _ = ax.set_xlabel("Median Intensity")
            _ = ax.set_xlim(xlim)
            _ = ax.set_ylim(ylim)
            _ = ax.legend()
            pl.tight_layout()

            return(fig, ax)

    def plot_volcano(self, ref_id, threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot", method="GMM_pvalue"):
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
        barplot: plot p-value bars instead of lines
        """
        # Extract fasta and positions
        start, end = self.__get_positions(ref_id)

        try:
            ref_pos_dict = self.results.query('ref_id==@ref_id').set_index('pos').to_dict('index')
        except NameError:
            raise NanocomporeError("It looks like there's not results slot in SampCompDB")

        # Make a list with all methods available
        methods=list(self.results)
        if method not in self._pvalue_tests:
            raise NanocomporeError("Method %s is not in the results dataframe. Please chose one of %s "%(method, self._pvalue_tests))

        # Parse line position per position
        d = OrderedDict()


        rp=self[ref_id]
        for pos in range(start, end+1):
            # Collect results for position
            res_dict = OrderedDict()
            if pos in ref_pos_dict:
                for k,v in ref_pos_dict[pos].items():
                    if k == method:
                        if not np.isnan(v):
                            res_dict[k] = -np.log10(v)
                        else:
                            res_dict[k] = 0
                res_dict['ref_kmer'] = ref_pos_dict[pos]['ref_kmer']
                if 'GMM_model' in rp[pos]['txComp']:
                    res_dict['LOR'] = rp[pos]['txComp']['GMM_model'][1]
                else:
                    res_dict['LOR'] = 0
            d[pos] = res_dict


        # Cast collected results to dataframe
        df = pd.DataFrame.from_dict(d, orient="index")
        if df.empty:
            raise NanocomporeError("No data available for the selected interval")
        # Define plotting style
        with pl.style.context(plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            _ = sns.scatterplot(x=df.LOR, y=df[method], palette=palette, ax=ax)
            for line in df.index.values:
                if df[method][line] > -np.log10(threshold):

                    ax.text(df.LOR[line], df[method][line], line, horizontalalignment='left', size='medium', color='black')
            _ = ax.axhline(y=-np.log10(threshold), color="grey", linestyle=":", label="pvalue={}".format(threshold))
            _ = ax.legend()
            _ = ax.set_ylabel("-log (pvalue)")
            _ = ax.set_xlabel("Log Odds Ratio")
            _ = ax.set_title(ref_id)
            pl.tight_layout()
            return(fig, ax)
