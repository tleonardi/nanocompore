#Native packages
import os, collections, sys

#3rd party packages
import numpy as np
import pandas as pd
from loguru import logger
from bedparse import bedline
from statsmodels.stats.multitest import multipletests

#Local Packages
from nanocompore.common import *
import nanocompore.SampComp_SQLDB as SampCompDb


class resultsManager():
    def __init__ (self, config):
        self._outpath = config.get_outpath()
        self._prefix = config.get_outprefix()
        self._overwrite = config.get_overwrite()
        self._bed_fn = config.get_bed()
        self._correction_method = config.get_correction_method()
        
        self._db = SampCompDb.SampComp_DB(outpath=self._outpath, prefix=self._prefix, overwrite=self._overwrite)
    
    def saveData(self, transcript, test_results, table=''):
        self._db.store_test_results(tx_name=transcript, test_results=test_results[0], table=table)
        self._db.store_read_level_data(tx_name=transcript, read_data=test_results[1])

    def closeDB(self):
        self._db.closeDB()

    def finish(self, valid_transcripts=set(), bed=False, bedgraph=False, pvalue_threshold=0):
        tests, shift_stats = self._getStatsTests()
        logger.debug("Gathering all the data for reporting")
        data = self._db.getAllData()
        logger.debug("Data gathered")
        data = self._addGenomicPositions(data, valid_transcripts=valid_transcripts)
        logger.debug("Added genomic positions to data")
        data = self._correct_pvalues(data, method=self._correction_method)
        logger.debug("Corrected pvalues")
        self._writeResultsTSV(data, tests)
        logger.debug("Wrote results TSV output file")
        self._writeResultsShiftStats(data, shift_stats)
        logger.debug("Wrote shift stats output file")
        if (bed or bedgraph) and self._bed_fn:
            if bed:
                self._writeBed(data, valid_transcripts=valid_transcripts, pvalue_threshold=pvalue_threshold)
            if bedgraph:
                self._writeBedgraph(data, valid_transcripts=valid_transcripts, pvalue_theshold=pvalue_threshold)
        elif (bed or bedgraph) and not self._bed_fn:
            raise NanocomporeError('Writing a bed file requires an input bed to map transcriptome coordinates to genome coordinates')
        self.closeDB()

    def _writeResultsTSV(self, data, stats_tests):
        #writes the results of all the statistical tests to a tsv file
        out_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_results.tsv")
        logger.debug(f"Starting to write results to {out_tsv}")
        columns_to_save = ['pos', 'chr', 'genomicPos', 'name', 'strand', 'sequence'] + ','.join(stats_tests).split(',')
        aliases = ['pos', 'chr', 'genomicPos', 'ref_id', 'strand', 'ref_kmer'] + ','.join(stats_tests).split(',')
        data.to_csv(out_tsv, sep='\t', columns=columns_to_save, header=aliases, index=False)

    def _writeResultsShiftStats(self, data, shift_stats):
        #writes the shift stats to a tsv file
        shift_stats_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_shift_stats.tsv")
        logger.debug(f"Starting to write results to {shift_stats_tsv}")
        columns_to_save = ['name', 'pos'] + ','.join(shift_stats).split(',')
        aliases = ['ref_id', 'pos'] + ','.join(shift_stats).split(',')
        data.to_csv(shift_stats_tsv, sep='\t', columns=columns_to_save, header=aliases, index=False)

    def _addGenomicPositions(self, data, valid_transcripts, convert='', assembly=''):
        '''
        if convert not in [None, "ensembl_to_ucsc", "ucsc_to_ensembl"]:
            sys.stderr.write("Convert value not valid\n")
            #raise NanocomporeError("Convert value not valid")
        if convert is not None and assembly is None:
            sys.stderr.write("The assembly argument is required in order to do the conversion. Choose one of 'hg38' or 'mm10'\n")
            #raise NanocomporeError("The assembly argument is required in order to do the conversion. Choose one of 'hg38' or 'mm10' ")
        '''
        if self._bed_fn:
            bed_annot={}
            try:
                with open(self._bed_fn) as tsvfile:
                    for line in tsvfile:
                        record_name=line.split('\t')[3]
                        if(record_name in valid_transcripts):
                            bed_annot[record_name]=bedline(line.split('\t'))

            except:
                raise NanocomporeError("Can't open BED file")
            if len(bed_annot) != len(valid_transcripts):
                raise NanocomporeError("Some references are missing from the BED file provided")

            data['genomicPos'] = data.apply(lambda row: bed_annot[row['name']].tx2genome(coord=row['pos'], stranded=True),axis=1)
            # This is very inefficient. We should get chr and strand only once per transcript, ideally when writing the BED file
            data['chr'] = data.apply(lambda row: bed_annot[row['name']].chr,axis=1)
            #if convert == "ensembl_to_ucsc":
            #    data['chr'] = data.apply(lambda row: bed_annot[row['name']].translateChr(assembly=assembly, target="ucsc", patches=True),axis=1)
            #elif convert == "ucsc_to_ensembl":
            #    data['chr'] = data.apply(lambda row: bed_annot[row['name']].translateChr(assembly=assembly, target="ens", patches=True),axis=1)
            data['strand'] = data.apply(lambda row: bed_annot[row['name']].strand,axis=1)
            logger.debug("Added genomic positions to the data")

        else:
            data['genomicPos'] = 'NA'
            data['chr'] = 'NA'
            data['strand'] = 'NA'
            logger.debug("No bed file was provided, genomic positions are not used and 'NA' will be used instead")
        
        return data

    def _writeBed(self, data, pvalue_threshold=0, span=5, title="Nanocompore Significant Sites"):
        if span < 1:
            #sys.stderr.write("span has to be >=1\n")
            raise NanocomporeError("span has to be >=1")

        data['end_pos'] = data['pos'] + span
        for test in data.columns:
            if 'pvalue' in  test:
                out_bed = os.path.join(self._outpath, f"{self._prefix}sig_sites_{test}_thr_{pvalue_threshold}.bed")
                logger.debug(f"Starting to write results to {out_bed}\n")
                columns_to_save = ['chr', 'genomicPos', 'end_pos', 'name', test, 'strand']
                data[data[test] <= pvalue_threshold].to_csv(out_bed, sep='\t', columns=columns_to_save, header=False, index=False)

    def _writeBedgraph(self, data, pvalue_threshold=0, title="Nanocompore Significant Sites"):
        data[data['strand'] == '+']['genomicPos'] = data[data['strand'] == '+']['genomicPos'] + 2
        data[data['strand'] == '+']['end_pos'] = data[data['strand'] == '+']['genomicPos'] + 3
        data[data['strand'] == '-']['genomicPos'] = data[data['strand'] == '-']['genomicPos'] - 2
        data[data['strand'] == '-']['end_pos'] = data[data['strand'] == '-']['genomicPos'] - 1

        for test in data.columns:
            if 'pvalue' in test:
                out_bedgraph = os.path.join(self._outpath, f"{self._prefix}sig_sites_{test}_thr_{pvalue_threshold}.bed")
                logger.debug(f"Starting to write results to {out_bedgraph}")
                columns_to_save = ['chr', 'genomicPos', 'end_pos', test,]
                data[data[test] <= pvalue_threshold].to_csv(out_bedgraph, sep='\t', columns=columns_to_save, header=False, index=False)     

    def _sort_headers_list(self, headers):
        headers.sort(key=lambda x: (not 'pvalue' in x, x))
        return headers  

    def _sort_shift_stat_headers(self, strings):
        conditions = ['c1', 'c2']
        stats = ['mean', 'median', 'sd']
        values = ['intensity', 'dwell']
        target = []
        for v in values:
            for s in stats:
                for c in conditions:
                    target.append(f"{c}_{s}_{v}")

        target_set = set(target)
        strings = [s for s in strings if s in target_set]
        strings.sort(key=lambda x: target.index(x) if x in target else len(target))
        return strings

    def _getStatsTests(self):
        tests = []
        shift_stats = []
        headers = self._db.get_kmer_table_headers()
        for header in headers:
            if self._any_header_in_string(['pvalue', 'LOR', 'GMM', 'logit', 'cluster_counts'], header):
                tests.append(header)
            elif self._any_header_in_string(['c1', 'c2'], header):
                shift_stats.append(header)
        return self._sort_headers_list(tests), self._sort_shift_stat_headers(shift_stats)

    def _any_header_in_string(self, check_list, string):
        for substring in check_list:
            if substring in string:
                return True
        return False

    def _correct_pvalues(self, df, method='fdr_bh'):
        for column in df.columns:
            if 'pvalue' in column:
                # Get the p-values from the column
                pvals = df[column].values

                # Replace NaN values with a large number so that they are not considered for correction
                nan_indices = np.isnan(pvals)
                pvals[nan_indices] = 1.0

                # Correct the p-values using the Benjamini-Hochberg method
                corrected_pvals = multipletests(pvals, method=method)[1]

                # Replace the original p-values with the corrected p-values in the dataframe
                df[column] = corrected_pvals
                
                # Replace the NaN values back in the dataframe
                df.loc[nan_indices, column] = np.nan
                logger.debug(f"pvalues for {column} have been corrected using {method}")
        return df