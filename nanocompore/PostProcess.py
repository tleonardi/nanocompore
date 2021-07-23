# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from loguru import logger

# Third party
# ...

# Local package
from nanocompore.common import *
from nanocompore.DataStore import DataStore_EventAlign, DataStore_SampComp

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class PostProcess(object):
    """Helper class for post-processing `SampComp` results"""

    def __init__(self, sampcomp_db_path:str, eventalign_db_path:str, bed_path:str=None):
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

        with DataStore_SampComp(self._sampcomp_db_path) as sc_db, \
             DataStore_EventAlign(self._eventalign_db_path) as ea_db:
            # do we have GMM results?
            query = "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'gmm_stats'"
            sc_db.cursor.execute(query)
            with_gmm = sc_db.cursor.fetchone() is not None
            query = "SELECT * FROM kmer_stats LEFT JOIN transcripts ON transcriptid = transcripts.id"
            if with_gmm:
                query += " LEFT JOIN gmm_stats ON kmer_stats.id = gmm_stats.kmer_statsid"
            query += " ORDER BY transcriptid, kmer"
            first_row = True
            shift_stat_columns = []
            univariate_pvalue_columns = []
            gmm_pvalue_columns = []
            for row in sc_db.cursor.execute(query):
                # retrieve k-mer sequence:
                ea_query = "SELECT sequence FROM kmers LEFT JOIN reads ON readid = reads.id WHERE transcriptid = ? AND position = ? LIMIT 1"
                ea_db.cursor.execute(ea_query, (row["transcriptid"], row["kmer"]))
                seq = ea_db.cursor.fetchone()[0]
                out_dict = {"transcript": row["name"],
                            "position": row["kmer"],
                            "sequence": seq}
                # TODO: add chromosome, genomic pos., strand information (from where?)
                if first_row: # check which columns we have (do this only once)
                    univariate_pvalue_columns = [col for col in row.keys()
                                                 if ("intensity_pvalue" in col) or ("dwell_pvalue" in col)]
                    if include_shift_stats:
                        shift_stat_columns = [col for col in row.keys() if col.startswith(("c1_", "c2_"))]
                    if with_gmm:
                        gmm_pvalue_columns = [col for col in row.keys() if "test_pvalue" in col]

                for col in shift_stat_columns:
                    out_dict[col] = row[col]
                for col in univariate_pvalue_columns:
                    out_dict[col] = row[col]
                if with_gmm:
                    out_dict["GMM_n_components"] = row["n_components"]
                    out_dict["GMM_cluster_counts"] = row["cluster_counts"]
                    out_dict["GMM_test_stat"] = row["test_stat"]
                    for col in gmm_pvalue_columns:
                        out_dict[col.replace("test", "GMM", 1)] = row[col]

                if first_row: # write header line
                    fp.write("\t".join(out_dict.keys()) + "\n")
                # write output data:
                fp.write("\t".join(str(x) for x in out_dict.values()) + "\n")
                first_row = False
