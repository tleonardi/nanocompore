import os

import numpy as np
import pandas as pd
import polars as pl

from loguru import logger
from statsmodels.stats.multitest import multipletests
from gtfparse import read_gtf

from nanocompore.common import decode_kmer
from nanocompore.database import ResultsDB


class Postprocessor():
    def __init__ (self, config):
        self._outpath = config.get_outpath()
        self._prefix = config.get_outprefix()
        self._result_exists_strategy = config.get_result_exists_strategy()
        self._correction_method = config.get_correction_method()
        self._config = config
        self._db = ResultsDB(config, init_db=False)


    def __call__(self):
        logger.info("Creating database indices.")
        self._db.index_database()

        all_columns = self._db.get_result_column_names()

        # The test_result_columns are all columns from the database
        # that need to be added to the result TSV.
        test_result_columns = self._get_test_result_columns(all_columns)
        self._export_results_tsv(test_result_columns)

        if self._config.get_export_shift_stats():
            logger.info("Writing the shift stats to a file.")
            shift_stats_columns = self._get_shift_stats_columns()
            self._export_shift_stats(shift_stats_columns)
        else:
            logger.info("Exporting shift statistics to a TSV is disabled. "
                        "However, you can still find them in the result database.")


    def _export_results_tsv(self, test_result_columns, chunksize=1_000_000):
        out_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_results.tsv")
        kmer_len = self._config.get_kit().len

        map_genomic_positions = bool(self._config.get_gtf())
        if map_genomic_positions:
            transcript_positions = self._get_transcript_positions(self._config.get_gtf())
            genomic_positions = np.vectorize(transcript_positions.get)

        # In order to do the multiple test correction we need
        # all p-values. However, we'd like to avoid loading
        # all data in memory. Hence, we read only the p-value
        # columns from the database and then the rest of the
        # results are processed in batches. I.e. we read each
        # batch, add the corresponding p-values and q-values
        # and write the append the batch to the results TSV.

        # We read the p-value columns we in memory and get add
        # the corrected q-values to them. We also get the
        # auto_test column if auto test was used, because we
        # need it in order to do the multiple test correction
        # independently for the different tests.
        pval_columns = [c for c in test_result_columns if 'pvalue' in c]
        p_values = self._db.get_columns(pval_columns)
        auto_test = None
        if 'auto_pvalue' in test_result_columns:
            auto_test = self._db.get_columns(['auto_test'])['auto_test']
        pq_values = self._correct_pvalues(p_values,
                                          method=self._correction_method,
                                          auto_test=auto_test)

        # The remaining columns will be loaded in batches from the db and written
        # to the output TSV.
        non_pqval_cols = [c for c in test_result_columns if 'pvalue' not in c]

        renamings = {'name': 'ref_id', 'kmer': 'ref_kmer'}

        base_columns = ['pos', 'name', 'kmer']

        # Since we have the p-values in memory,
        # we only want to get the base columns and
        # the non-p-value test result columns.
        select_columns = base_columns + non_pqval_cols
        test_columns = list(pq_values.columns) + non_pqval_cols

        output_columns = ['pos', 'chr', 'genomicPos', 'ref_id', 'strand', 'ref_kmer'] + test_columns

        logger.info(f"Writing the test results to: {out_tsv}")

        completed_rows = 0
        for chunk in self._db.get_results(select_columns, chunksize=chunksize):
            nrows = chunk.shape[0]

            chunk['genomicPos'] = 'NA'
            chunk['chr'] = 'NA'
            chunk['strand'] = 'NA'

            # Make sure counts are integer
            chunk = chunk.astype({col: 'int'
                                  for col in chunk.columns
                                  if col.endswith('_mod') or col.endswith('_unmod')})

            chunk['kmer'] = [decode_kmer(code, kmer_len) for code in chunk['kmer']]

            chunk.rename(columns=renamings, inplace=True)
            qvals_chunk = pq_values[completed_rows:completed_rows+nrows].reset_index(drop=True)
            chunk = pd.concat([chunk, qvals_chunk], axis=1)

            if map_genomic_positions:
                genomic_data = genomic_positions(chunk.ref_id.str.split('|').str[0])
                coords_getter = np.vectorize(
                        lambda tx, pos: (tx['seqname'],
                                         tx['strand'],
                                         tx['exons'](pos)))
                chrom, strand, gpos = coords_getter(genomic_data, chunk.pos)
                chunk['chr'] = chrom
                chunk['strand'] = strand
                chunk['genomicPos'] = gpos

            mode = 'w' if completed_rows == 0 else 'a'
            header = completed_rows == 0
            chunk.to_csv(out_tsv, sep='\t', columns=output_columns, index=False, header=header, mode=mode)

            completed_rows += nrows


    def _export_shift_stats(self, shift_stats_columns, chunksize=1_000_000):
        shift_stats_tsv = os.path.join(self._outpath, f"{self._prefix}nanocompore_shift_stats.tsv")
        base_columns = ['name', 'pos']
        renamings = {'name': 'ref_id'}
        select_columns = base_columns + shift_stats_columns
        for i, chunk in enumerate(self._db.get_results(select_columns, chunksize=chunksize)):
            chunk.rename(columns=renamings, inplace=True)
            chunk = chunk.round(decimals=4)
            mode = 'w' if i == 0 else 'a'
            header = i == 0
            chunk.to_csv(shift_stats_tsv, sep='\t', index=False, header=header, mode=mode)


    def _get_transcript_positions(self, gtf_path):
        gtf = read_gtf(gtf_path, features=['exon'])
        gtf = gtf.select('transcript_id',
                         'seqname',
                         'strand',
                         pl.struct('start', 'end').alias('exons'))\
                 .group_by('transcript_id', 'seqname', 'strand')\
                 .agg('exons')\
                 .to_pandas()\
                 .set_index('transcript_id')\
                 .to_dict('index')
        for transcript in gtf.values():
            transcript['exons'] = self._get_genomic_mapper(transcript['exons'],
                                                           transcript['strand'])
        return gtf


    def _get_genomic_mapper(self, exons, strand):
        if strand == '+':
            positions = np.concatenate([np.arange(e['start'], e['end']+1) for e in exons])
        else:
            positions = np.concatenate([np.arange(e['end'], e['start']-1, -1) for e in exons])
        # GTF files use 1-based indexing.
        # Hence, we subtract one to convert to 0-based indices.
        positions -= 1
        return np.vectorize(dict(zip(range(len(positions)),
                                     positions)).get)


    def _sort_headers_list(self, headers):
        return sorted(headers,
                      key=lambda x: (not ('pvalue' in x or 'qvalue' in x),
                                     ('_mod' in x or '_unmod' in x),
                                     x))


    def _get_test_result_columns(self, columns):
        keywords = ['pvalue', 'qvalue', 'LOR', 'GMM', 'logit', 'auto_test', '_mod', '_unmod']
        selected =  [col
                     for col in columns
                     if any(kw in col for kw in keywords)]
        return self._sort_headers_list(selected)


    def _get_shift_stats_columns(self):
        conditions = ['c1', 'c2']
        stats = ['mean', 'median', 'sd']
        values = ['intensity', 'dwell']
        return [f"{c}_{s}_{v}"
                for v in values
                for s in stats
                for c in conditions]


    def _correct_pvalues(self, df, method='fdr_bh', auto_test=None):
        for column in df.columns:
            if column == 'auto_pvalue' and method == 'fdr_bh':
                pvals = np.array(df[column].values, dtype=float)
                all_qvals = np.full(pvals.shape, np.nan, dtype=float)

                logger.debug(f"Starting to correct pvalues for {column} with {method}")
                # We need to correct the p-values for a total of N
                # tests. However, those p-values were obtained by
                # performing different tests and are not directly
                # comparable. Using the standard Benjamini-Hochberg
                # will introduce bias against tests that are less
                # sensitive by pushing their p-values to the bottom
                # of the ranking and then applying more severe penalties
                # for them. As a simple workaround, we perform the
                # correction separately (thus ordering different
                # tests independently) and then we multiply by
                # a correction coefficient to compensate for the
                # fact that the correction was done for less than
                # N tests.
                for test_type in auto_test.unique():
                    test_mask = auto_test == test_type
                    # We'll ignore NaNs when performing the correction.
                    qvals = self._multipletests_filter_nan(pvals[test_mask], method)
                    if len(qvals) > 0:
                        qvals *= len(pvals) / len(qvals)
                    all_qvals[test_mask] = qvals
                df['auto_qvalue'] = all_qvals
                logger.debug(f"pvalues for {column} have been corrected using {method}")
            elif 'pvalue' in column:
                pvals = np.array(df[column].values, dtype=float)

                # We'll ignore NaNs when performing the correction.
                logger.debug(f"Starting to correct pvalues for {column} with {method}")
                qvals = self._multipletests_filter_nan(pvals, method)

                # Replace the original p-values with the corrected p-values in the dataframe
                corrected_column = column.replace("pvalue", "qvalue")
                df[corrected_column] = qvals

                logger.debug(f"pvalues for {column} have been corrected using {method}")

        df.sort_index(axis=1, inplace=True)
        return df


    @staticmethod
    def _multipletests_filter_nan(pvalues: np.ndarray, method="fdr_bh") -> np.ndarray:
        """
        Performs p-value correction for multiple hypothesis testing
        using the method specified. The pvalues list can contain
        np.nan values, which are ignored during p-value correction.

        Parameters
        ----------
        pvalues : np.ndarray
            Array of p-values to correct.
        method :
            Which correction strategy to use. Should be supported
            by statsmodels's multipletests method.

        Returns
        -------
        np.ndarray
            Array of q-values.

        Example
        -------
        input = np.array([0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan])
        Postprocessor._multipletests_filter_nan(input)
        array([0.13333333, 0.016, nan, 0.016, 0.5, 0.45714286, 0.016, 0.008, nan, nan, 0.016, nan])
        """
        nans = np.isnan(pvalues)
        corrected_p_values = multipletests(pvalues[~nans], method=method)[1]
        qvalues = np.copy(pvalues)
        qvalues[~nans] = corrected_p_values
        return qvalues

