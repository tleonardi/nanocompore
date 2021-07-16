# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from loguru import logger

# Third party
from statsmodels.stats.multitest import multipletests


#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class PostProcess(object):
    """Helper class for post-processing `SampComp` results"""

    def __init___(self, sampcomp_db_path:str, eventalign_db_path:str, bed_path:str=None):
        self._sampcomp_db_path = sampcomp_db_path
        self._eventalign_db_path = eventalign_db_path
        self._bed_path = bed_path


    def save_all(self, outpath_prefix=None, pvalue_thr=0.01):
        """
        Save all text reports including genomic coordinate if a bed file was provided
        * outpath_prefix
            outpath + prefix to use as a basename for output files.
            If not given, it will use the same prefix as the database.
        * pvalue_thr
            pvalue threshold to report significant sites in bed files
        """
        if not outpath_prefix:
            outpath_prefix = self._db_path.replace("SampComp.db", "")
        logger.debug("Save reports to {}".format(outpath_prefix))

        # Save reports
        logger.debug("\tSaving extended tabular report")
        self.save_report(output_fn = outpath_prefix + "nanocompore_results.tsv")
        logger.debug("\tSaving shift results")
        self.save_shift_stats(output_fn = outpath_prefix + "nanocompore_shift_stats.tsv")

        # Save bed and bedgraph files for each method used
        if self._bed_path:
            logger.debug("\tSaving significant genomic coordinates in Bed and Bedgraph format")
            for m in self._metadata["pvalue_tests"]:
                self.save_to_bed(
                    output_fn = outpath_prefix+"sig_sites_{}_thr_{}.bed".format(m, pvalue_thr),
                    bedgraph=False, pvalue_field=m, pvalue_thr=pvalue_thr, span=5, title="Nanocompore Significant Sites")
                self.save_to_bed(
                    output_fn = outpath_prefix+"sig_sites_{}_thr_{}.bedgraph".format(m, pvalue_thr),
                    bedgraph=True, pvalue_field=m, pvalue_thr=pvalue_thr, title="Nanocompore Significant Sites")


    def save_to_bed(self, output_fn=None, bedgraph=False, pvalue_field=None, pvalue_thr=0.01, span=5, convert=None, assembly=None, title=None):
        """
        Save the position of significant positions in the genome space in BED6 or BEDGRAPH format.
        The resulting file can be used in a genome browser to visualise significant genomic locations.
        The option is only available if `SampCompDB` if initialised with a BED file containing genome annotations.
        * output_fn
            Path to file where to write the data
        * bedgraph
            save file in bedgraph format instead of bed
        * pvalue_field
            specifies what column to use as BED score (field 5, as -log10)
        * pvalue_thr
            only report positions with pvalue<=thr
        * span
            The size of each BED feature.
            If size=5 (default) features correspond to kmers.
            If size=1 features correspond to the first base of each kmer.
        * convert
            one of 'ensembl_to_ucsc' or 'ucsc_to_ensembl". Convert chromosome named between Ensembl and Ucsc conventions
        * assembly
            required if convert is used. One of "hg38" or "mm10"
        """
        if self._bed_path is None:
            raise NanocomporeError("In order to generate a BED file PostProcess needs to be initialised with a transcriptome BED")
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
                    bed_file.write('track type=bed name="%s" description="%s"\n' % (title, pvalue_field))
                else:
                    bed_file.write('track type=bedGraph name="%s" description="%s"\n' % (title, pvalue_field))

            Record = namedtuple('Record', ['chr', 'genomicPos', 'ref_id', 'strand', 'ref_kmer', pvalue_field])
            threshold = -log(pvalue_thr, 10)
            for record in self.results[list(Record._fields)].itertuples(index=False, name="Record"):
                pvalue = getattr(record, pvalue_field)
                if np.isnan(pvalue):
                    pvalue = 0
                elif pvalue < sys.float_info.min:
                    pvalue = -log(sys.float_info.min, 10)
                else:
                    pvalue = -log(pvalue, 10)
                if not bedgraph and pvalue < threshold:
                    continue
                if bedgraph:
                    if record.strand == "+":
                        start_pos = record.genomicPos + 2
                    else:
                        start_pos = record.genomicPos - 2
                    end_pos = start_pos + 1
                else:
                    if record.strand == "+":
                        start_pos = record.genomicPos
                    else:
                        start_pos = record.genomicPos - span + 1
                    end_pos = start_pos + span
                line = bedline([record.chr, start_pos, end_pos, "%s_%s" % (record.ref_id, record.ref_kmer),
                                pvalue, record.strand])
                if convert == "ensembl_to_ucsc":
                    line = line.translateChr(assembly=assembly, target="ucsc", patches=True)
                elif convert == "ucsc_to_ensembl":
                    line = line.translateChr(assembly=assembly, target="ens", patches=True)
                if bedgraph:
                    bed_file.write("%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end, line.score))
                else:
                    bed_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (line.chr, line.start, line.end,
                                                                 line.name, line.score, line.strand))


    def save_report(self, output_fn:str=None, include_shift_stats:bool=True):
        """
        Saves a tabulated text dump of the database containing all the statistical results for all the positions
        * output_fn
            Path to file where to write the data. If None, data is returned to the standard output.
        """
        ## TODO: can this be done in a "with ..." clause?
        if output_fn is None:
            fp = sys.stdout
        elif isinstance(output_fn, str):
            try:
                fp = open(output_fn, "w")
            except:
                raise NanocomporeError("Error opening output file %s" % output_fn)
        else:
            raise NanocomporeError("output_fn needs to be a string or None")

        shift_stat_columns = []
        if include_shift_stats:
            shift_stat_columns = ["c1_mean_intensity", "c2_mean_intensity",
                                  "c1_median_intensity", "c2_median_intensity",
                                  "c1_sd_intensity", "c2_sd_intensity",
                                  "c1_mean_dwell", "c2_mean_dwell",
                                  "c1_median_dwell", "c2_median_dwell",
                                  "c1_sd_dwell", "c2_sd_dwell"]

        with DataStore_SampComp(self._sampcomp_db_path) as sc_db, \
             DataStore_EventAlign(self._eventalign_db_path) as ea_db:
            # Which statistical tests were performed?
            query = "SELECT DISTINCT test FROM univariate_results"
            univar_tests = [row["test"] for row in sc_db.cursor.execute(query)]
            query = "SELECT DISTINCT test FROM gmm_results"
            gmm_tests = [row["test"] for row in sc_db.cursor.execute(query)]
            # Generate headers
            headers = ['pos', 'chr', 'genomicPos', 'ref_id', 'strand', 'ref_kmer']
            for test in sorted(univar_tests):
                headers += [f"{test}_dwell_pvalue", f"{test}_intensity_pvalue"]
            if gmm_tests:
                # TODO: what if GMM was fitted, but no test were performed?
                headers += ["GMM_cov_type", "GMM_n_clust", "cluster_counts"]
                if "logit" in gmm_tests:
                    headers += ["GMM_logit_pvalue", "Logit_LOR"]
                if "anova" in gmm_tests:
                    headers += ["GMM_anova_pvalue", "Anova_delta_logit"]
            # Write headers to file
            fp.write('\t'.join([str(i) for i in headers]) + '\n')

            # Merge kmer information with transcript name:
            columns = ["kmer_stats.id", "transcriptid", "kmer AS pos", "name AS ref_id"] + shift_stat_columns
            columns = ", ".join(columns)
            query = f"SELECT {columns} FROM kmer_stats LEFT JOIN transcripts ON transcriptid = transcripts.id ORDER BY transcriptid, kmer"
            for row in sc_db.cursor.execute(query):
                db_data = dict(row)
                # Get p-values etc.:
                id = db_data["id"]
                if univar_tests:
                    query = f"SELECT test, intensity_pvalue, dwell_pvalue FROM univariate_results WHERE kmer_statsid = {id}"
                    for row2 in sc_db.cursor.execute(query):
                        test = row2["test"]
                        db_data[test + "_intensity_pvalue"] = row2["intensity_pvalue"]
                        db_data[test + "_dwell_pvalue"] = row2["dwell_pvalue"]
                if gmm_tests:
                    query = f"SELECT test, test_pvalue, test_stat FROM gmm_results WHERE gmm_statsid = {id}"
                    for row2 in sc_db.cursor.execute(query):
                        test = row2["test"]
                        db_data[test + "_intensity_pvalue"] = row2["intensity_pvalue"]
                        db_data[test + "_dwell_pvalue"] = row2["dwell_pvalue"]



            # TODO: where does chromosome and genomic pos. information come from?


        # We loop over the IDs so that ref_pos_list can be prefetched for each transcript
        for cur_id in self.ref_id_list:
            cur_ref_pos_list = self[cur_id]
            for record in self.results[self.results.ref_id == cur_id ].itertuples():
                if "GMM" in self._metadata["comparison_methods"]:
                    record_txComp = cur_ref_pos_list[record.pos]['txComp']
                line = []
                for f in headers:
                    if f in record._fields:
                        line.append(getattr(record, f))
                    elif f == "GMM_cov_type":
                        line.append(record_txComp['GMM_model']['model'].covariance_type)
                    elif f == "GMM_n_clust":
                        line.append(record_txComp['GMM_model']['model'].n_components)
                    elif f == "cluster_counts":
                        line.append(record_txComp['GMM_model']['cluster_counts'])
                    elif f == "Anova_delta_logit":
                        line.append(record_txComp['GMM_anova_model']['delta_logit'])
                    elif f == "Logit_LOR":
                        line.append(record_txComp['GMM_logit_model']['coef'])
                    else: line.append("NA")
                fp.write('\t'.join([ str(i) for i in line ])+'\n')
        fp.close()


    def save_shift_stats(self, output_fn=None):
        """
        Save the mean, median and sd intensity and dwell time for each condition and for each position.
        This can be used to evaluate the intensity of the shift for significant positions.
        * output_fn
            Path to file where to write the data. If None, data is returned to the standard output.
        """
        if output_fn is None:
            fp = sys.stdout
        elif isinstance(output_fn, str):
            try:
                fp = open(output_fn, "w")
            except:
                raise NanocomporeError("Error opening output file %s" % output_fn)
        else:
            raise NanocomporeError("output_fn needs to be a string or None")

        headers = ['c1_mean_intensity', 'c2_mean_intensity', 'c1_median_intensity', 'c2_median_intensity', 'c1_sd_intensity', 'c2_sd_intensity', 'c1_mean_dwell', 'c2_mean_dwell', 'c1_median_dwell', 'c2_median_dwell', 'c1_sd_dwell', 'c2_sd_dwell']
        fp.write('\t'.join([ str(i) for i in ["ref_id", "pos"]+headers ])+'\n')
        for tx, refpos in self:
            for pos, refpos_list in enumerate(refpos):
                if "txComp" in refpos_list:
                    ss = refpos_list['txComp']['shift_stats']
                    if list(ss.keys()) != headers:
                        raise NanocomporeError("Mismatch in shift_stats headers")
                    line = [tx, pos, *ss.values()]
                    fp.write('\t'.join([ str(i) for i in line ])+'\n')
        fp.close()


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
        if all([np.isnan(p) for p in pvalues]):
            return pvalues

        pvalues_no_nan = [p for p in pvalues if not np.isnan(p)]
        corrected_p_values = multipletests(pvalues_no_nan, method=method)[1]
        for i, p in enumerate(pvalues):
            if np.isnan(p):
                corrected_p_values = np.insert(corrected_p_values, i, np.nan, axis=0)
        return(corrected_p_values)
