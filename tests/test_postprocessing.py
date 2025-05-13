import copy
import os
import tempfile
import shutil

import numpy as np
import pandas as pd

from nanocompore.postprocessing import Postprocessor
from nanocompore.config import Config

from tests.common import BASIC_CONFIG


def test_multipletests_filter_nan():
    pvals = np.array([0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan])
    expected_qvals = np.array([0.13333333, 0.016, np.nan, 0.016, 0.5, 0.45714286, 0.016, 0.008, np.nan, np.nan, 0.016, np.nan])

    qvals = Postprocessor._multipletests_filter_nan(pvals)

    assert np.allclose(qvals, expected_qvals, equal_nan=True)


def test_export_results_tsv_single_chunk_no_genomic_positions():
    config_yaml = copy.deepcopy(BASIC_CONFIG)

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_results.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        all_columns = postprocessor._db.get_result_column_names()
        test_result_columns = postprocessor._get_test_result_columns(all_columns)
        postprocessor._export_results_tsv(test_result_columns)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 21)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['pos',
                         'chr',
                         'genomicPos',
                         'ref_id',
                         'strand',
                         'ref_kmer',
                         'GMM_chi2_pvalue',
                         'GMM_chi2_qvalue',
                         'KS_dwell_pvalue',
                         'KS_dwell_qvalue',
                         'KS_intensity_pvalue',
                         'KS_intensity_qvalue',
                         'GMM_LOR',
                         'KD1_mod',
                         'KD1_unmod',
                         'KD2_mod',
                         'KD2_unmod',
                         'WT1_mod',
                         'WT1_unmod',
                         'WT2_mod',
                         'WT2_unmod']
        assert np.all(results.columns == expected_cols)
        assert np.all(np.isnan(results.genomicPos))
        assert np.all(np.isnan(results.chr))
        assert np.all(np.isnan(results.strand))


def test_export_results_tsv_single_chunk_with_genomic_positions():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['gtf'] = 'tests/fixtures/test_annotation.gtf'

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_results.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        all_columns = postprocessor._db.get_result_column_names()
        test_result_columns = postprocessor._get_test_result_columns(all_columns)
        postprocessor._export_results_tsv(test_result_columns)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 21)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['pos',
                         'chr',
                         'genomicPos',
                         'ref_id',
                         'strand',
                         'ref_kmer',
                         'GMM_chi2_pvalue',
                         'GMM_chi2_qvalue',
                         'KS_dwell_pvalue',
                         'KS_dwell_qvalue',
                         'KS_intensity_pvalue',
                         'KS_intensity_qvalue',
                         'GMM_LOR',
                         'KD1_mod',
                         'KD1_unmod',
                         'KD2_mod',
                         'KD2_unmod',
                         'WT1_mod',
                         'WT1_unmod',
                         'WT2_mod',
                         'WT2_unmod']
        assert np.all(results.columns == expected_cols)

        test_position = results[(results.pos == 89) & (results.ref_id == 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')].iloc[0]
        assert test_position.genomicPos == 5529654
        assert test_position.chr == 'chr7'
        assert test_position.strand == '-'

        test_position = results[(results.pos == 1982) & (results.ref_id == 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|')].iloc[0]
        assert test_position.genomicPos == 5527184
        assert test_position.chr == 'chr7'
        assert test_position.strand == '-'


def test_export_results_tsv_multiple_chunks_no_genomic_positions():
    config_yaml = copy.deepcopy(BASIC_CONFIG)

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_results.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        all_columns = postprocessor._db.get_result_column_names()
        test_result_columns = postprocessor._get_test_result_columns(all_columns)
        postprocessor._export_results_tsv(test_result_columns, chunksize=100)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 21)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['pos',
                         'chr',
                         'genomicPos',
                         'ref_id',
                         'strand',
                         'ref_kmer',
                         'GMM_chi2_pvalue',
                         'GMM_chi2_qvalue',
                         'KS_dwell_pvalue',
                         'KS_dwell_qvalue',
                         'KS_intensity_pvalue',
                         'KS_intensity_qvalue',
                         'GMM_LOR',
                         'KD1_mod',
                         'KD1_unmod',
                         'KD2_mod',
                         'KD2_unmod',
                         'WT1_mod',
                         'WT1_unmod',
                         'WT2_mod',
                         'WT2_unmod']
        assert np.all(results.columns == expected_cols)
        assert np.all(np.isnan(results.genomicPos))
        assert np.all(np.isnan(results.chr))
        assert np.all(np.isnan(results.strand))


def test_export_results_tsv_multple_chunks_with_genomic_positions():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['gtf'] = 'tests/fixtures/test_annotation.gtf'

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_results.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        all_columns = postprocessor._db.get_result_column_names()
        test_result_columns = postprocessor._get_test_result_columns(all_columns)
        postprocessor._export_results_tsv(test_result_columns, chunksize=100)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 21)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['pos',
                         'chr',
                         'genomicPos',
                         'ref_id',
                         'strand',
                         'ref_kmer',
                         'GMM_chi2_pvalue',
                         'GMM_chi2_qvalue',
                         'KS_dwell_pvalue',
                         'KS_dwell_qvalue',
                         'KS_intensity_pvalue',
                         'KS_intensity_qvalue',
                         'GMM_LOR',
                         'KD1_mod',
                         'KD1_unmod',
                         'KD2_mod',
                         'KD2_unmod',
                         'WT1_mod',
                         'WT1_unmod',
                         'WT2_mod',
                         'WT2_unmod']
        assert np.all(results.columns == expected_cols)

        test_position = results[(results.pos == 89) & (results.ref_id == 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')].iloc[0]
        assert test_position.genomicPos == 5529654
        assert test_position.chr == 'chr7'
        assert test_position.strand == '-'

        test_position = results[(results.pos == 1982) & (results.ref_id == 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|')].iloc[0]
        assert test_position.genomicPos == 5527184
        assert test_position.chr == 'chr7'
        assert test_position.strand == '-'


def test_export_shift_stats_single_chunk():
    config_yaml = copy.deepcopy(BASIC_CONFIG)

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_shift_stats.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        shift_columns = postprocessor._get_shift_stats_columns()
        postprocessor._export_shift_stats(shift_columns)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 14)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['ref_id',
                         'pos',
                         'c1_mean_intensity',
                         'c2_mean_intensity',
                         'c1_median_intensity',
                         'c2_median_intensity',
                         'c1_sd_intensity',
                         'c2_sd_intensity',
                         'c1_mean_dwell',
                         'c2_mean_dwell',
                         'c1_median_dwell',
                         'c2_median_dwell',
                         'c1_sd_dwell',
                         'c2_sd_dwell']
        assert np.all(results.columns == expected_cols)

        test_position = results[(results.pos == 89) & (results.ref_id == 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')].iloc[0]
        assert test_position.c1_mean_intensity == 12116.667
        assert test_position.c1_mean_dwell == 1.7454

        test_position = results[(results.pos == 1982) & (results.ref_id == 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|')].iloc[0]
        assert test_position.c1_mean_intensity == 7186.5
        assert test_position.c1_mean_dwell == 1.3476


def test_export_shift_stats_multiple_chunks():
    config_yaml = copy.deepcopy(BASIC_CONFIG)

    # Since we may run tests in parallel we want to
    # avoid having two tests write at the same path
    # simultaneously.
    with tempfile.TemporaryDirectory() as tmp:
        config_yaml['outpath'] = tmp

        shutil.copyfile('tests/fixtures/out_sampComp_sql.db',
                        os.path.join(tmp, 'out_sampComp_sql.db'))

        results_path = os.path.join(tmp, 'out_nanocompore_shift_stats.tsv')

        config = Config(config_yaml)
        postprocessor = Postprocessor(config)

        shift_columns = postprocessor._get_shift_stats_columns()
        postprocessor._export_shift_stats(shift_columns, chunksize=100)

        results = pd.read_csv(results_path, sep='\t')

        assert results.shape == (3437, 14)
        assert set(results.ref_id.unique()) == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                                                'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
        expected_cols = ['ref_id',
                         'pos',
                         'c1_mean_intensity',
                         'c2_mean_intensity',
                         'c1_median_intensity',
                         'c2_median_intensity',
                         'c1_sd_intensity',
                         'c2_sd_intensity',
                         'c1_mean_dwell',
                         'c2_mean_dwell',
                         'c1_median_dwell',
                         'c2_median_dwell',
                         'c1_sd_dwell',
                         'c2_sd_dwell']
        assert np.all(results.columns == expected_cols)

        test_position = results[(results.pos == 89) & (results.ref_id == 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')].iloc[0]
        assert test_position.c1_mean_intensity == 12116.667
        assert test_position.c1_mean_dwell == 1.7454

        test_position = results[(results.pos == 1982) & (results.ref_id == 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|')].iloc[0]
        assert test_position.c1_mean_intensity == 7186.5
        assert test_position.c1_mean_dwell == 1.3476

