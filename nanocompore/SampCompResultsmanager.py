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
import nanocompore.SampComp_SQLDB as SampCompDB


class resultsManager():
    def __init__ (self, config):
        self._outpath = config.get_outpath()
        self._prefix = config.get_outprefix()
        self._result_exists_strategy = config.get_result_exists_strategy()
        self._bed_fn = config.get_bed()
        self._correction_method = config.get_correction_method()
        self._config = config

        self._db = SampCompDB.SampCompDB(config, init_db=True)


    def save_results(self, transcript, test_results):
        self._db.save_test_results(transcript, test_results)


    def finish(self, valid_transcripts=set(), bed=False, bedgraph=False, pvalue_threshold=0):
        self._db.index_database()
        logger.debug("Created all table indexes")
        logger.debug("Gathering all the data for reporting")
        data = self._db.get_all_results()
        tests, shift_stats = self._getStatsTests(data)
        logger.debug("Data gathered")
        data = self._addGenomicPositions(data, valid_transcripts=valid_transcripts)
        data['kmer'] = data['kmer'].apply(lambda encoding: decode_kmer(encoding, self._config.get_kit().len))
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


    def filter_already_processed_transcripts(self, transcripts):
        existing_transcripts = set(self._db.get_transcripts())
        return [tx for tx in transcripts if tx.ref_id not in existing_transcripts]


    def _writeResultsTSV(self, data, stats_tests):
        #writes the results of all the statistical tests to a tsv file
        out_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_results.tsv")
        logger.debug(f"Starting to write results to {out_tsv}")
        columns_to_save = ['pos', 'chr', 'genomicPos', 'name', 'strand', 'kmer'] + ','.join(stats_tests).split(',')
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

    def _getStatsTests(self, data):
        tests = []
        shift_stats = []
        for header in data.columns:
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
                pvals = np.array(df[column].values, dtype=float)

                # Correct the p-values using the Benjamini-Hochberg method.
                # We'll ignore NaNs when performing the correction.
                logger.debug(f"Starting to correct pvalues for {column} with {method}")
                corrected_pvals = self.__multipletests_filter_nan(pvals, method)

                # Replace the original p-values with the corrected p-values in the dataframe
                df[column] = corrected_pvals

                logger.debug(f"pvalues for {column} have been corrected using {method}")
        return df


    @staticmethod
    def __multipletests_filter_nan(pvalues, method="fdr_bh"):
        """
        Performs p-value correction for multiple hypothesis testing
        using the method specified. The pvalues list can contain
        np.nan values, which are ignored during p-value correction.
        test: input=[0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan]
        out: array([0.13333333, 0.016, nan, 0.016, 0.5, 0.45714286, 0.016, 0.008, nan, nan, 0.016, nan])
        """
        if np.isnan(pvalues).all():
            return pvalues

        pvalues_no_nan = [p for p in pvalues if not np.isnan(p)]
        corrected_p_values = multipletests(pvalues_no_nan, method=method)[1]
        current_corrected = 0
        results = np.empty(len(pvalues))
        for i, p in enumerate(pvalues):
            if np.isnan(p):
                results[i] = np.nan
            else:
                results[i] = corrected_p_values[current_corrected]
                current_corrected += 1
        return results

