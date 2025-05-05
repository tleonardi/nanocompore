import os

import numpy as np
import pandas as pd

from bedparse import bedline
from loguru import logger
from statsmodels.stats.multitest import multipletests

from nanocompore.common import NanocomporeError
from nanocompore.common import decode_kmer
from nanocompore.database import ResultsDB


class Postprocessor():
    def __init__ (self, config):
        self._outpath = config.get_outpath()
        self._prefix = config.get_outprefix()
        self._result_exists_strategy = config.get_result_exists_strategy()
        self._bed_fn = config.get_bed()
        self._correction_method = config.get_correction_method()
        self._config = config
        self._db = ResultsDB(config, init_db=False)


    def __call__(self):
        logger.info("Creating database indices.")
        self._db.index_database()

        logger.info("Gathering all the data for reporting.")
        data = self._db.get_all_results()

        # Make sure counts are integer
        data = data.astype({col: 'int'
                            for col in data.columns
                            if col.endswith('_mod') or col.endswith('_unmod')})

        logger.info("Adding genomic positions to data.")
        data = self._add_genomic_positions(data)

        logger.info("Decoding kmers.")
        kmer_len = self._config.get_kit().len
        data['kmer'] = [decode_kmer(code, kmer_len) for code in data['kmer']]

        logger.info(f'Correcting p-values using "{self._correction_method}"')
        data = self._correct_pvalues(data, method=self._correction_method)

        logger.info("Writing the test results to a file.")
        tests, shift_stats = self._get_stats_tests(data)
        self._write_results_tsv(data, tests)

        logger.info("Writing the shift stats to a file.")
        self._write_results_shift_stats(data, shift_stats)

        if self._bed_fn:
            data.sort_values(['chr', 'strand', 'name', 'genomicPos'], inplace=True)
            logger.info("Writing the results to a bed file.")
            self._write_bed(data)


    def _write_results_tsv(self, data, stats_tests):
        #writes the results of all the statistical tests to a tsv file
        out_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_results.tsv")
        logger.debug(f"Starting to write results to {out_tsv}")
        columns_to_save = ['pos', 'chr', 'genomicPos', 'name', 'strand', 'kmer']
        aliases = ['pos', 'chr', 'genomicPos', 'ref_id', 'strand', 'ref_kmer']
        columns_to_save.extend(stats_tests)
        aliases.extend(stats_tests)
        data.to_csv(out_tsv, sep='\t', columns=columns_to_save, header=aliases, index=False)


    def _write_results_shift_stats(self, data, shift_stats):
        """
        Writes the shift statistics to a tsv file.
        """
        shift_stats_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_shift_stats.tsv")
        logger.debug(f"Starting to write results to {shift_stats_tsv}")
        columns_to_save = ['name', 'pos']
        aliases = ['ref_id', 'pos']
        columns_to_save.extend(shift_stats)
        aliases.extend(shift_stats)
        data.to_csv(shift_stats_tsv, sep='\t', columns=columns_to_save, header=aliases, index=False)


    def _add_genomic_positions(self, data, convert='', assembly=''):
        if self._bed_fn:
            bed_annot = {}
            try:
                with open(self._bed_fn) as tsvfile:
                    for line in tsvfile:
                        record_name = line.split('\t')[3]
                        bed_annot[record_name] = bedline(line.split('\t'))

            except Exception:
                raise NanocomporeError("Can't open BED file")


            data['genomicPos'] = pd.Series(dtype='int')
            data['chr'] = pd.Series(dtype='str')
            data['strand'] = pd.Series(dtype='str')
            for i, row in data.iterrows():
                transcript = row['name'].split('|')[0]
                annotation = bed_annot[transcript]
                genomic_pos = annotation.tx2genome(row['pos'], stranded=True)
                data.loc[i, 'genomicPos'] = int(genomic_pos)
                data.loc[i, 'chr'] = annotation.chr
                data.loc[i, 'strand'] = annotation.strand
            logger.debug("Added genomic positions to the data")
        else:
            data['genomicPos'] = 'NA'
            data['chr'] = 'NA'
            data['strand'] = 'NA'
            logger.debug("No bed file was provided, genomic positions are not used and 'NA' will be used instead")

        return data


    def _write_bed(self, data):
        kit = self._config.get_kit()
        # The interval calculation is a bit tricky, so here's an example:
        # Suppose we use RNA002 with kmer = 5 and center position (most infulential base)
        # at the 4th base.
        # If genomicPos is 3 (that's the fourth base, because it's a 0-based index).
        # chromosome positions: 0 1 2 3 4 5
        #                             ^
        #                    this is the kmer center
        # We want to obtain the interval [0, 5), because BED files use
        # left-inclusive, right-exclusive intervals.
        # So the proper calculation would be:
        # start_pos = genomicPos - center + 1 = 3 - 4  + 1 = 0
        # end_pos = genomicPos + kmer_len - center + 1 = 3 + 5 - 4 + 1 = 5
        data['start_pos'] = data['genomicPos'] - kit.center + 1
        data['end_pos'] = data['genomicPos'] + kit.len - kit.center + 1
        for test in data.columns:
            if 'qvalue' in test:
                out_bed = os.path.join(self._outpath, f"{self._prefix}sig_sites_{test}.bed")
                logger.debug(f"Starting to write results to {out_bed}\n")
                columns_to_save = ['chr', 'start_pos', 'end_pos', 'name', test, 'strand']
                data.to_csv(out_bed,
                            sep='\t',
                            columns=columns_to_save,
                            header=False,
                            index=False)


    def _sort_headers_list(self, headers):
        return sorted(headers,
                      key=lambda x: (not ('pvalue' in x or 'qvalue' in x),
                                     ('_mod' in x or '_unmod' in x),
                                     x))


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


    def _get_stats_tests(self, data):
        tests = []
        shift_stats = []
        for header in data.columns:
            keywords = ['pvalue', 'qvalue', 'LOR', 'GMM', 'logit', 'auto_test', '_mod', '_unmod']
            if any([kw in header for kw in keywords]):
                tests.append(header)
            elif any([kw in header for kw in ['c1', 'c2']]):
                shift_stats.append(header)
        return self._sort_headers_list(tests), self._sort_shift_stat_headers(shift_stats)


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
                corrected_column = column.replace("pvalue", "qvalue")
                df[corrected_column] = corrected_pvals

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

