# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import OrderedDict, namedtuple
import shelve
from dbm import error as dbm_error
from math import log
import re
import sys

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
from nanocompore.common import counter_to_str, access_file, NanocomporeError
from nanocompore import models

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class SampCompDB (object):
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
            with shelve.open (db_fn, flag='r') as db:
                # Try to get metadata from db
                try:
                    metadata = db['__metadata']
                    self._comparison_method = metadata['comparison_method']
                    self._sequence_context = metadata['sequence_context']
                    self._min_coverage = metadata['min_coverage']
                    self._n_samples = metadata['n_samples']
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
            self._model_dict = models.RNA_model_dict
        else:
            raise NanocomporeError ("Only RNA is implemented at the moment")

        self.bed_fn = bed_fn

        # Create results DF with adjusted p-values
        self.__calculate_results()

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


    #~~~~~~~~~~~~~~PRIVATE  METHODS~~~~~~~~~~~~~~#

    def __calculate_results(self, adjust=True, methods=None):
        # Compose a lists with the name of the results
        tests=[]
        if methods is None:
            methods = self._comparison_method
        for m in methods:
            if m in ["MW", "KS", "TT"]:
                tests.append("{}_intensity_pvalue".format(m))
                tests.append("{}_dwell_pvalue".format(m))
                if self._sequence_context:
                    tests.append("{}_intensity_pvalue_context_{}".format(m, self._sequence_context))
                    tests.append("{}_dwell_pvalue_context_{}".format(m, self._sequence_context))
            elif m == "GMM":
                tests.append("GMM_pvalue")
                if self._sequence_context:
                    tests.append("GMM_pvalue_context_{}".format(self._sequence_context))
        tests.append("lowCov")
        # We open the DB rather that calling __getitem__ to avoid
        # opening and closing the file for every transcript
        with shelve.open(self._db_fn, flag = "r") as db:
            df = pd.DataFrame([dict({x:y for a,b in v.items() if a == "txComp" for x,y in b.items()  if x in tests}, ref=ref_id, pos=k, ref_kmer=v['ref_kmer'])  for ref_id, rec in db.items() for k,v in rec.items() if ref_id!="__metadata" ])

        if self.bed_fn:
            bed_annot={}
            try:
                with open(self.bed_fn) as tsvfile:
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
            df=df[['ref', 'pos', 'chr', 'strand', 'genomicPos', 'ref_kmer']+tests]
        else:
            df=df[['ref', 'pos', 'ref_kmer']+tests]

        if adjust:
            for col in tests:
                if "pvalue" in col:
                    df['adjusted_'+col] = self.__multipletests_filter_nan(df[col], method="fdr_bh")
        self.results = df


    def __get_fasta_seq (self, ref_id, start, end):
        """ Extract fasta record corresponding to ref with error handling """
        try:
            with Fasta(self._fasta_fn) as fasta:
                fasta = (fasta [ref_id])
                return str (fasta[start:end])
        except KeyError:
            raise NanocomporeError("Reference id not present in fasta file")


    def __get_positions (self, ref_id, start=None, end=None):
        """ Verify start and end and if not available try to infer them"""
        try:
            with Fasta(self._fasta_fn) as fasta:
                max_len = len(fasta [ref_id])
        except KeyError:
            raise NanocomporeError("Reference id not present in fasta file")

        if not start:
            start = 0
        if not end:
            end = max_len
        if start > end:
            raise NanocomporeError ("End coordinate has to be higher or equal to start")
        if start < 0:
            raise NanocomporeError ("Coordinates have to be higher that 0")
        if end > max_len:
            raise NanocomporeError ("Coordinates have to be lower than the ref_id sequence length ({})".format(max_len))

        return (start, end)


    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def save_to_bed (self, output_fn, bedgraph=False, pvalue_field=None, pvalue_thr=0.01, span=5, convert=None, assembly=None, title=None):
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

            Record = namedtuple('Record', ['chr', 'genomicPos', 'ref','strand', 'ref_kmer', pvalue_field ])
            for record in self.results[ list(Record._fields) ].itertuples(index=False, name="Record"):
                pvalue = getattr(record, pvalue_field)
                if np.isnan(pvalue): 
                    pvalue=0
                else:
                    pvalue=-log(pvalue, 10)
                if not bedgraph and pvalue >= -log(pvalue_thr, 10):
                    line=bedline([record.chr, record.genomicPos, record.genomicPos+span, f"{record.ref}_{record.ref_kmer}", pvalue, record.strand])
                    if convert is "ensembl_to_ucsc":
                        line=line.translateChr(assembly=assembly, target="ucsc", patches=True)
                    elif convert is "ucsc_to_ensembl":
                        line=line.translateChr(assembly=assembly, target="ens", patches=True)
                    bed_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.name, line.score, line.strand))
                elif bedgraph:
                    line=bedline([record.chr, record.genomicPos+2, record.genomicPos+3, f"{record.ref}_{record.ref_kmer}", pvalue, record.strand])
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
        r = re.compile("adjusted_*")
        methods = list(filter(r.match, list(self.results)))
        
        headers = ['chr', 'pos', 'ref','strand', 'ref_kmer', 'lowCov']+methods
        # Read extra GMM info from the shelve
        if "GMM" in self._comparison_method:
            headers += ['LOR', 'clusters']
            gmm_info=OrderedDict()
            for tx, refpos in self:
                gmm_info[tx] = {k:{'lor': v['txComp']['GMM_model'][1], 'clusters':v['txComp']['GMM_model'][3]} for k,v in refpos.items() if "GMM_model" in v['txComp']}
        fp.write('\t'.join([ str(i) for i in headers ])+'\n')
        for record in self.results[['chr', 'pos', 'ref','strand', 'ref_kmer', 'lowCov']+methods].values.tolist():
            if "GMM" in self._comparison_method:
                try:
                    lor = gmm_info[record[2]][record[1]]['lor']
                    clusters = gmm_info[record[2]][record[1]]['clusters']
                    clusters = '#'.join([ ','.join([str(x) for x in i]) for i in clusters])
                    record += [lor, clusters]
                except KeyError:
                    record += ["nan", "nan"]
            fp.write('\t'.join([ str(i) for i in record ])+'\n')
        fp.close()
    


    def list_most_significant_positions (self, n=10):
        pass

    def list_most_significant_references (self, n=10):
        pass

    #~~~~~~~~~~~~~~PLOTTING METHODS~~~~~~~~~~~~~~#
    def plot_pvalue(self, ref_id, start=None, end=None, threshold=0.01, figsize=(30,10), palette="Set2", plot_style="ggplot", method=None, barplot=False):
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
        start, end = self._SampCompDB__get_positions(ref_id, start, end)

        try:
            ref_pos_dict = self.results.query('ref==@ref_id').set_index('pos').to_dict('index')
        except NameError:
            raise NanocomporeError("It looks like there's not results slot in SampCompDB")

        # Make a list with all methods available
        methods=list(self.results)

        # If method not provided, set it to a default regex
        if method is None:
            method="adjusted_*"
        # Parse the method as regex if string
        if isinstance(method, str):
            r = re.compile(method)
            method = list(filter(r.match, methods))
        elif not isinstance(method, list):
            raise NanocomporeError("Method must be either a string or a list")

        for m in method:
            if m not in methods:
                raise NanocomporeError("Method %s is not in the results dataframe"%m)

        # Parse line position per position
        d = OrderedDict()

        for pos in range(start, end+1):
            # Collect results for position
            res_dict = OrderedDict ()
            if pos in ref_pos_dict:
                for k,v in ref_pos_dict[pos].items():
                    if k in method:
                        if not np.isnan(v):
                            res_dict[k] = -np.log10(v)
                        else:
                            res_dict[k] = 0
                res_dict['ref_kmer'] = ref_pos_dict[pos]['ref_kmer']
            d[pos] = res_dict


        # Cast collected results to dataframe
        df = pd.DataFrame.from_dict(d, orient="index")
        if df.empty:
            raise NanocomporeError("No data available for the selected interval")
        # Define plotting style
        with pl.style.context(plot_style):
            fig, ax = pl.subplots(figsize=figsize)
            if not barplot:
                _ = sns.lineplot(data=df.drop('ref_kmer', axis=1), palette=palette, ax=ax, dashes=False)
            else:
                df = df.reset_index().melt(id_vars=["index","ref_kmer"], var_name="method", value_name="pvalue").set_index('index')
                _ = sns.barplot(x=df.index.values, y=df.pvalue, hue=df.method,  palette=palette, ax=ax)
            _ = ax.axhline(y=-np.log10(threshold), color="grey", linestyle=":", label="pvalue={}".format(threshold))
            _ = ax.legend()
            _ = ax.set_ylabel("-log (pvalue)")
            _ = ax.set_xlabel("Reference position")
            _ = ax.set_title(ref_id)
            if not barplot:
                if end-start<30:
                    _ = ax.set_xticks(df.index.values)
                    _ = ax.set_xticklabels( [ "{}\n{}".format(i[0], i[1]) for i in zip(df.index, df.ref_kmer) ] )
            else:
                if end-start<30:
                    _ = ax.set_xticklabels( [ "{}\n{}".format(i[0], i[1]) for i in zip(df.index, df.ref_kmer) ] )
                else:
                    step = (end-start)//10
                    breaks = list(range(start, end+1, step))
                    _ = ax.set_xticklabels( [ i if i in breaks else "" for i in df.index ]  )
            pl.tight_layout()
            return(fig, ax)


    def plot_signal (self, ref_id, start=None, end=None, split_samples=False, feature="intensity", figsize=(30,10), palette="Set2", plot_style="ggplot", bw=0.25):
        """
        Plot the dwell time and median intensity distribution position per positon in a split violin plot representation.
        It is pointless to plot more than 50 positions at once as it becomes hard to distiguish
        ref_id: Valid reference id name in the database
        start: Start coordinate (Must be higher or equal to 0)
        end: End coordinate (included) (must be lower or equal to the reference length)
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
        start, end = self.__get_positions (ref_id, start, end)

        # Parse line position per position
        l_intensity = []
        l_dwell = []
        valid=0
        x_ticks_list = []
        model_means_list = []

        # Extract data from database if position in db
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                valid+=1
                ref_kmer=ref_data[pos]['ref_kmer']
                x_ticks_list.append("{}\n{}".format(pos, ref_kmer))
                model_means_list.append(self._model_dict[ref_kmer][0])

                for k1, v1 in ref_data[pos]['data'].items():
                    # Collect dwell or median data
                    for k2, v2 in v1.items():
                        lab = "{}_{}".format(k1, k2) if split_samples else k1

                        for value in v2["intensity"]:
                            l_intensity.append ((pos, lab, value))
                        for value in v2["dwell"]:
                            l_dwell.append ((pos, lab, np.log10(value)))

            else:
                l_intensity.append ((pos, None, None))
                l_dwell.append ((pos, None, None))
                x_ticks_list.append(str(pos))
                model_means_list.append(None)

        # Check that we found valid position and cast collected results to dataframe
        if not valid:
            raise NanocomporeError ("No data available for selected coordinates")
        df_intensity = pd.DataFrame (l_intensity, columns=["pos", "lab", "value"])
        df_dwell = pd.DataFrame (l_dwell, columns=["pos", "lab", "value"])

        # Define ploting style
        with pl.style.context (plot_style):
            fig, (ax1, ax2) = pl.subplots(2,1, figsize=figsize, sharex=True)

            # Plot median intensity violin + model mean
            _ = sns.violinplot (x="pos", y="value", hue="lab", data=df_intensity, ax=ax1, split=not split_samples, inner="quartile", bw=bw, linewidth=1, scale="area", palette=palette)
            _ = ax1.plot (model_means_list, color="black", marker="x", label="Model Mean", linestyle="")
            _ = ax1.set_ylabel ("Mean Intensity")
            _ = ax1.set_xlabel ("")
            _ = ax1.legend ()

            _ = sns.violinplot ( x="pos", y="value", hue="lab", data=df_dwell, ax=ax2, split=not split_samples, inner="quartile", bw=bw, linewidth=1, scale="area", palette=palette)
            _ = ax2.set_ylabel ("log10 (Dwell Time)")
            _ = ax2.set_xlabel ("Reference position")
            _ = ax2.set_xlim (-1, end-start+1)
            _ = ax2.set_xticklabels (x_ticks_list)
            _ = ax2.legend ()

            # Adjust display
            _ = fig.suptitle("Reference:{}  Start:{}  End:{}".format(ref_id, start, end), y=1.01)

            pl.tight_layout()
            return (fig, (ax1, ax2))

    def plot_coverage (self, ref_id, start=None, end=None, scale=False, figsize=(30,5), palette="Set2", plot_style="ggplot"):
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
        """
        # Extract data for ref_id
        ref_data = self[ref_id]
        start, end = self.__get_positions (ref_id, start, end)

        # Parse line position per position
        l=[]
        valid=0
        # Extract data from database if position in db
        for pos in np.arange (start, end+1):
            if pos in ref_data:
                valid+=1
                for k1, v1 in ref_data[pos]['data'].items():
                    for k2, v2 in v1.items():
                        l.append ((pos, "{}_{}".format(k1, k2), v2["coverage"]))
            else:
                l.append ((pos, None, None))

        # Check that we found valid position and cast collected results to dataframe
        if not valid:
            raise NanocomporeError ("No data available for selected coordinates")
        df = pd.DataFrame (l, columns=["pos", "Sample", "cov"])
        if scale:
            df['cov'] = df.groupby('Sample')['cov'].apply(lambda x: x/max(x))

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
            if not scale:
                _ = ax.axhline  (y=self._min_coverage, linestyle=":", color="grey", label="minimal coverage")
            _ = ax.set_ylim (0, None)
            _ = ax.set_xlim (start, end)
            _ = ax.set_title ("Reference:{}  Start:{}  End:{}".format(ref_id, start, end))
            _ = ax.set_ylabel ("Coverage")
            _ = ax.set_xlabel ("Reference position")
            _ = ax.legend()

            pl.tight_layout()
            return (fig, ax)

    def plot_position(self, ref_id, pos=None, split_samples=False, figsize=(30,10), palette="Set2",  plot_style="ggplot", xlim=None, ylim=None, alpha=0.3, pointSize=20, scatter=True, kde=True, model=False, gmm_levels=50):
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
        if pos not in ref_data:
            raise NanocomporeError("No data available for the selected position")

        # Extract data from database if position in db
        ref_kmer = ref_data[pos]['ref_kmer']
        data = ref_data[pos]['data']

        # Disretize palette
        if split_samples:
            colors = sns.mpl_palette(palette, self._n_samples)
        else:
            colors = sns.mpl_palette(palette, 2)

        # Collect and transform data in dict
        plot_data_dict = OrderedDict ()
        i = 0
        for cond_label, cond_data in data.items():
            if split_samples:
                for rep_lab, rep_data in cond_data.items():
                    plot_data_dict["{}_{}".format(cond_label, rep_lab)] = {
                        "intensity":scale(rep_data["intensity"]),
                        "dwell":scale(np.log10(rep_data["dwell"])),
                        "color":colors[i]}
                    i+=1
            else:
                intensity_list = []
                dwell_list = []
                for rep_lab, rep_data in cond_data.items():
                    intensity_list.append(rep_data["intensity"])
                    dwell_list.append(rep_data["dwell"])
                plot_data_dict[cond_label] = {
                    "intensity":scale(np.concatenate(intensity_list)),
                    "dwell":scale(np.log10(np.concatenate(dwell_list))),
                    "color":colors[i]}
                i+=1

        if model:
            model = self[ref_id][pos]['txComp']['GMM_model'][4]
            if isinstance(model, GaussianMixture):
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
        with pl.style.context(plot_style):
            # Plot dwell and median
            fig, ax = pl.subplots(figsize=figsize)

            for label, d in plot_data_dict.items():
                if kde:
                    _ = sns.kdeplot(
                        d["intensity"],
                        d["dwell"],
                        cmap=sns.light_palette(d["color"], as_cmap=True),
                        ax=ax,
                        clip=((min(d["intensity"]), max(d["intensity"])), (min(d["dwell"]),max(d["dwell"]))))
                if scatter:
                    _ = ax.scatter (
                        x=d["intensity"],
                        y=d["dwell"],
                        color=d["color"],
                        label=label,
                        alpha=alpha,
                        s=pointSize)
            if model:
                _ = ax.contour(X, Y, Z, levels=gmm_levels, alpha=alpha, colors="black")
            # Adjust display
            _ = ax.set_title ("%s\n%s (%s)"%(ref_id,pos, ref_kmer))
            _ = ax.set_ylabel ("log10 (Dwell Time)")
            _ = ax.set_xlabel ("Median Intensity")
            if xlim:
                _ = ax.set_xlim(xlim)
            if ylim:
                _ = ax.set_ylim(ylim)
            _ = ax.legend()
            pl.tight_layout()

            return (fig, ax)

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
