# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
import os.path

# Third party
from loguru import logger

# Local package
from nanocompore.common import *
from nanocompore.DataStore import DataStore_master, DataStore_transcript

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
class PostProcess(object):
    """Helper class for post-processing `SampComp` results"""

    def __init__(self, input_dir:str, master_db:str="nanocompore.db", bed_path:str=None):
        self._input_dir = input_dir
        self._master_db_path = os.path.join(input_dir, master_db)
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


    def save_report(self, output_fn:str=None, significance_thresholds={"adj_gmm_pvalue": 0.01}, details=False):
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

        where = []
        order = []
        for k, v in significance_thresholds.items():
            where.append(f"{k} <= {v}")
            order.append(k)
        if details or not order:
            order = ["transcriptid", "kmer_pos"]
        # TODO: depending on 'details', results will be ordered by p-value or by transcript/pos. - problem?
        sql = "SELECT * FROM test_results LEFT JOIN (SELECT id, name, subdir FROM transcripts) ON transcriptid = id"
        if where:
            sql += f" WHERE " + " OR ".join(where)
        sql += " ORDER BY " + ", ".join(order)
        exclude_cols = ["transcriptid", "id", "name", "subdir"]
        with DataStore_master(self._master_db_path) as master:
            if not details: # write output from master DB directly to TSV
                master.cursor.execute(sql)
                row = master.cursor.fetchone()
                if row:
                    self.write_tsv_header(fp, row)
                    self.write_tsv_row(fp, row)
                    for row in master.cursor:
                        self.write_tsv_row(fp, row)
            else:
                # include GMM stats?
                master.cursor.execute("SELECT value FROM parameters WHERE step = 'SC' AND name = 'fit_gmm'")
                row = master.cursor.fetchone()
                with_gmm = (row is not None) and (row[0] == "True")
                master.cursor.execute(sql)
                row = master.cursor.fetchone()
                first_row = True # do we need to write the header?
                while row:
                    current_tx = row["name"]
                    db_path = os.path.join(self._input_dir, str(row["subdir"]), current_tx + ".db")
                    with DataStore_transcript(db_path, current_tx, row["transcriptid"]) as db:
                        while row and (row["name"] == current_tx): # still on the same transcript
                            # get data from transcript DB:
                            kmer_pos = row["kmer_pos"]
                            seq_query = "SELECT sequenceid FROM kmers WHERE position = ? LIMIT 1"
                            seq_row = db.cursor.execute(seq_query, (kmer_pos, )).fetchone()
                            seq = list(master.sequence_mapping.keys())[seq_row[0]]
                            kmer_query = "SELECT * FROM kmer_stats WHERE kmer_pos = ?"
                            kmer_row = db.cursor.execute(kmer_query, (kmer_pos, )).fetchone()
                            if with_gmm:
                                gmm_query = "SELECT cluster_counts FROM gmm_stats WHERE kmer_pos = ?"
                                gmm_row = db.cursor.execute(gmm_query, (kmer_pos, )).fetchone()
                                if not gmm_row: # may be 'None' if best GMM only had one component
                                    gmm_row = [None]
                            else:
                                gmm_row = None
                            if first_row:
                                self.write_tsv_header(fp, row, seq, kmer_row, gmm_row)
                                first_row = False
                            self.write_tsv_row(fp, row, seq, kmer_row, gmm_row)
                            row = master.cursor.fetchone()


    @staticmethod
    def write_tsv_header(fp, master_row, seq=None, kmer_row=None, gmm_row=None):
        fp.write("transcript\tkmer_pos")
        if seq: # squeeze kmer sequence in at this point
            fp.write("\tkmer_seq")
        for k in master_row.keys():
            if "pvalue" in k:
                fp.write("\t" + k)
        if kmer_row:
            for k in kmer_row.keys()[1:]: # skip 'kmer_pos'
                if "pvalue" not in k: # skip p-value columns (already written from master DB)
                    fp.write("\t" + k)
            if gmm_row:
                fp.write("\tcluster_counts")
        fp.write("\n")


    @staticmethod
    def write_tsv_row(fp, master_row, seq=None, kmer_row=None, gmm_row=None):
        fp.write(master_row["name"] + "\t" + str(master_row["kmer_pos"]))
        if seq: # squeeze kmer sequence in at this point
            fp.write("\t" + seq)
        for k in master_row.keys():
            if "pvalue" in k:
                fp.write("\t" + str(master_row[k]))
        if kmer_row:
            for k in kmer_row.keys()[1:]: # skip 'kmer_pos'
                if "pvalue" not in k: # skip p-value columns (already written from master DB)
                    fp.write("\t" + str(kmer_row[k]))
            if gmm_row:
                fp.write("\t" + str(gmm_row[0]))
        fp.write("\n")
