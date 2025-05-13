import copy

from unittest.mock import Mock, call

import torch
import numpy as np
import pandas as pd

from nanocompore.comparisons import TranscriptComparator
from nanocompore.comparisons import calculate_lors
from nanocompore.comparisons import get_contigency_matrices
from nanocompore.config import Config
from nanocompore.common import MOTOR_DWELL_POS

from tests.common import BASIC_CONFIG
from tests.common import MockWorker
from tests.common import naneq


def test_add_shift_stats():
    config = Config(BASIC_CONFIG)
    data = torch.tensor([
        # pos 0
        [
            # The columns are:
            # intensity, dwell
            [1.3, 1.1], # read 0
            [1.2, 0.1], # read 1
            [1.3, 0.3], # read 2
            [1.0, 0.2], # read 3
        ],
        # pos 1
        [
            [1.3, 1.1], # read 4
            [1.2, 0.1], # read 5
            [1.1, 0.2], # read 6
            [0.6, 0.4], # read 6
        ],
    ])
    conditions = torch.tensor([0, 1, 0, 1])


    comparator = TranscriptComparator(config, MockWorker())

    results = {}
    comparator._add_shift_stats(results,
                                data,
                                conditions,
                                'cpu')

    assert get_float(results['c1_mean_intensity'][0]) == 1.3
    assert get_float(results['c1_mean_intensity'][1]) == 1.2
    assert get_float(results['c1_mean_dwell'][0]) == 0.7
    assert get_float(results['c1_mean_dwell'][1]) == 0.65

    assert get_float(results['c2_mean_intensity'][0]) == 1.1
    assert get_float(results['c2_mean_intensity'][1]) == 0.9
    assert get_float(results['c2_mean_dwell'][0]) == 0.15
    assert get_float(results['c2_mean_dwell'][1]) == 0.25

    assert get_float(results['c1_sd_intensity'][1]) == 0.1


def test_add_shift_stats_with_motor_dwell():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['motor_dwell_offset'] = 11
    config = Config(config_yaml)
    data = torch.tensor([
        # pos 0
        [
            # The columns are:
            # intensity, dwell, motor_dwell
            [1.3, 1.1, 0.6], # read 0
            [1.2, 0.1, 1.0], # read 1
            [1.3, 0.3, 0.4], # read 2
            [1.0, 0.2, 0.3], # read 3
        ],
        # pos 1
        [
            [1.3, 1.1, 0.6], # read 4
            [1.2, 0.1, 1.0], # read 5
            [1.1, 0.2, 0.4], # read 6
            [0.6, 0.4, 0.2], # read 6
        ],
    ])
    conditions = torch.tensor([0, 1, 0, 1])

    comparator = TranscriptComparator(config, MockWorker())

    results = {}
    comparator._add_shift_stats(results,
                                data,
                                conditions,
                                'cpu')

    assert get_float(results['c1_mean_intensity'][0]) == 1.3
    assert get_float(results['c1_mean_intensity'][1]) == 1.2
    assert get_float(results['c1_mean_dwell'][0]) == 0.7
    assert get_float(results['c1_mean_dwell'][1]) == 0.65
    assert get_float(results['c1_mean_motor_dwell'][0]) == 0.5
    assert get_float(results['c1_mean_motor_dwell'][1]) == 0.5

    assert get_float(results['c2_mean_intensity'][0]) == 1.1
    assert get_float(results['c2_mean_intensity'][1]) == 0.9
    assert get_float(results['c2_mean_dwell'][0]) == 0.15
    assert get_float(results['c2_mean_dwell'][1]) == 0.25
    assert get_float(results['c2_mean_motor_dwell'][0]) == 0.65
    assert get_float(results['c2_mean_motor_dwell'][1]) == 0.6

    assert get_float(results['c1_sd_intensity'][1]) == 0.1


def test_nonparametric_test_KS():
    config = Config(BASIC_CONFIG)
    comparator = TranscriptComparator(config, MockWorker())
    data = torch.tensor([
        # pos 3
        [
            [97.5, np.log10(0.08), np.log10(0.3)], # read 1
            [98.8, np.log10(0.012), np.log10(0.03)], # read 2
            [84.3, np.log10(0.4), np.nan], # read 3
            [113.2, np.log10(0.06), np.nan], # read 4
            [np.nan, np.nan, np.nan], # read 5
            [np.nan, np.nan, np.nan] # read 6
            ],
        # pos 5
        [
            [98.5, np.log10(0.09), np.nan], # read 1
            [99.8, np.log10(0.013), np.nan], # read 2
            [85.3, np.log10(0.5), np.nan], # read 3
            [114.2, np.log10(0.07), np.nan], # read 4
            [np.nan, np.nan, np.nan], # read 5
            [np.nan, np.nan, np.nan] # read 6
            ],
        # pos 14 - motor kmer of pos 3
        [
            [82.3, np.log10(0.3), np.nan], # read 1
            [111.2, np.log10(0.03), np.nan], # read 2
            [np.nan, np.nan, np.nan], # read 3
            [np.nan, np.nan, np.nan], # read 4
            [96.5, np.log10(0.07), np.nan], # read 5
            [97.8, np.log10(0.01), np.nan] # read 6
            ]
        ], dtype=float)
    conditions = torch.tensor([0, 1, 0, 1, 0, 1])

    results = comparator._nonparametric_test('KS', data, conditions)

    assert 'KS_intensity_pvalue' in results
    assert 'KS_dwell_pvalue' in results

    assert len(results['KS_intensity_pvalue']) == 3
    assert len(results['KS_dwell_pvalue']) == 3


def test_combine_context_pvalues():
    yaml = copy.deepcopy(BASIC_CONFIG)
    yaml['sequence_context_weights'] = 'harmonic'
    yaml['sequence_context'] = 3
    config = Config(yaml)
    comparator = TranscriptComparator(config, MockWorker())
    results = pd.DataFrame({'pos': range(1, 17),
                            'GMM_chi2_pvalue': [0.001,
                                                np.nan,
                                                0.01,
                                                0.103,
                                                0.1,
                                                0.301,
                                                0.901,
                                                np.nan,
                                                0.001,
                                                0.701,
                                                np.nan,
                                                0.98,
                                                0.99,
                                                0.03,
                                                0.99,
                                                0.99]})
    # Remove the first two nan to simulate low coverage
    # positions that were not tested. The third nan is
    # left in place to simulate a tested position that
    # has failed.
    results.drop([1, 7], axis='index', inplace=True)

    actual = comparator._combine_context_pvalues(results, 'GMM_chi2_pvalue')
    actual = np.array(actual)

    # These values are the results taken from the
    # v1 implementation of the context merging.
    expected = np.array([0.002123656975661474,
                         0.0019018653049881877,
                         0.00711410035724805,
                         0.06601202462390741,
                         0.04439107441454254,
                         0.17963717597756168,
                         0.011383510423954581,
                         0.27956447545269386,
                         np.nan,
                         0.471918623078623,
                         0.8109385225797305,
                         0.3482236872871895,
                         0.8342987523763911,
                         0.9504729672694502])
    assert np.equal(np.nan_to_num(actual.round(3), -1),
                    np.nan_to_num(expected.round(3), -1)).all()


def test_calculate_lor():
    contingencies = torch.tensor([
        # pos 0
        [[13, 11],
         [12, 10]],
        # pos 1
        [[9, 11],
         [7, 21]],
        # pos 2
        [[2, 1],
         [7, 1]]])
    assert torch.all(calculate_lors(contingencies) == torch.tensor([-0.015, 0.898, -1.253]))


def test_get_contigency_matrices():
    # conditions has shape (reads,)
    conditions = torch.tensor([0, 1, 0, 1, 0, 1])
    # predictions has shape (positions, reads)
    predictions = torch.tensor([[0, 0, 0, 1, 1, 1],
                                [0, 1, 0, 1, 0, 1]])
    expected = torch.tensor([
        # pos 0 contengency matrix
        [[2.0, 1],
         [1, 2]],
        # pos 1 contingency matrix
        [[3, 0],
         [0, 3]]])
    assert torch.all(get_contigency_matrices(conditions, predictions) == expected)


def test_get_contigency_matrices_with_gaps():
    # conditions has shape (reads,)
    conditions = torch.tensor([0, 1, 0, 1, 0, 1])
    # predictions has shape (positions, reads)
    predictions = torch.tensor([[0, 0, 0, np.nan, 1, 1],
                                [0, 1, 0, np.nan, 0, 1]])
    expected = torch.tensor([
        # pos 0 contengency matrix
        [[2.0, 1],
         [1, 1]],
        # pos 1 contingency matrix
        [[3, 0],
         [0, 2]]])
    assert torch.all(get_contigency_matrices(conditions, predictions) == expected)


def test_get_cluster_counts():
    config = Config(BASIC_CONFIG)
    # Condition 0 is the knockdown (depleted condition).
    # E.g. for the first position:
    # Most of the points for the non-depleted condition (1)
    # are found in cluster 1, then we can guess that the
    # modification cluster is 1.
    # For the third, cluster 0 will be the modified one.
    conditions = torch.tensor([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
    samples = torch.tensor([0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3])
    # We have predictions for three positions.
    predictions = torch.tensor([[0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1],
                                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]])

    contingency = get_contigency_matrices(conditions, predictions)

    comparator = TranscriptComparator(config, MockWorker())

    cluster_counts = comparator._get_cluster_counts(contingency, samples, conditions, predictions)

    assert np.all(cluster_counts['kd1_mod'] == np.array([0, 0, 0]))
    assert np.all(cluster_counts['kd2_mod'] == np.array([1, 0, 0]))
    assert np.all(cluster_counts['wt1_mod'] == np.array([1, 2, 2]))
    assert np.all(cluster_counts['wt2_mod'] == np.array([3, 4, 4]))
    assert np.all(cluster_counts['kd1_unmod'] == np.array([2, 2, 2]))
    assert np.all(cluster_counts['kd2_unmod'] == np.array([2, 3, 3]))
    assert np.all(cluster_counts['wt1_unmod'] == np.array([1, 0, 0]))
    assert np.all(cluster_counts['wt2_unmod'] == np.array([1, 0, 0]))


def test_gmm_test():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['motor_dwell_offset'] = 1
    config = Config(config_yaml)
    comparator = TranscriptComparator(config, MockWorker())

    data = torch.tensor([
        # pos 0 (2d test because high proportion of missing motor dwell)
        [[1., 1., np.nan],
         [np.nan, np.nan, 1.0],
         [1., 1., np.nan],
         [np.nan, np.nan, 1.0],
         [1., 1., 1.0],
         [1., 1., np.nan]],
        # pos 1 (should be run with 3d test)
        [[2., 2., 2.0],
         [2., 2., 2.0],
         [2., 2., 2.0],
         [2.0, 2.0, 2.0],
         [1., 1., 2.0],
         [2., 2., 2.0]],
        # pos 2 (2d test because high proportion of missing motor dwell)
        [[np.nan, np.nan, np.nan],
         [np.nan, np.nan, 3.0],
         [3., 3., np.nan],
         [np.nan, np.nan, 3.0],
         [3., 3., np.nan],
         [np.nan, np.nan, np.nan]],
    ])
    samples = torch.arange(6)
    conditions = torch.tensor([0, 0, 0, 1, 1, 1])

    gmm_results_2d = {'GMM_chi2_pvalue': np.array([0.1, 0.82]),
                      'GMM_LOR': np.array([0.8, 0.2])}
    gmm_results_3d = {'GMM_chi2_pvalue': np.array([0.01]),
                      'GMM_LOR': np.array([1.3])}

    def test(data, samples, conditions, device):
        if data.shape[2] == 3:
            return gmm_results_3d
        return gmm_results_2d

    comparator._gmm_test_split = Mock(side_effect=test)

    results = comparator._gmm_test(data, samples, conditions, 'cpu')

    callargs = comparator._gmm_test_split.call_args_list

    # Validate that the gmm_test_split is called once for pos 1 with 3d data
    # and once for positions 0 and2 with 2d data. Since we don't want to
    # impose the order we try both ways.
    assert ((naneq(callargs[0].args[0], data[[True, False, True], :, :MOTOR_DWELL_POS])
             and
             naneq(callargs[1].args[0], data[[False, True, False]]))
            or
            (naneq(callargs[0].args[0], data[[False, True, False]])
             and
             naneq(callargs[1].args[0], data[[True, False, True], :, :MOTOR_DWELL_POS])))

    # Validate that the results of the 2d and 3d tests
    # are properly merged in the correct order.
    assert np.array_equal(results['GMM_chi2_pvalue'], np.array([0.1, 0.01, 0.82]))
    assert np.array_equal(results['GMM_LOR'], np.array([0.8, 1.3, 0.2]))


def test_gmm_test_split_single_component():
    config = Config(BASIC_CONFIG)
    comparator = TranscriptComparator(config, MockWorker())

    rand_gen = np.random.default_rng(seed=42)
    data = torch.tensor(rand_gen.standard_normal((1, 100, 2)))
    samples = torch.tensor([samp
                            for samp in [0, 1, 2, 3]
                            for _ in range(25)])
    conditions = torch.tensor([cond
                               for cond in [0, 1]
                               for _ in range(50)])

    # The data is sampled from a normal distribution
    # and one component should fit it better (i.e.
    # should have lower BIC value). We should've
    # detected that and omit reporting the p-value
    # from the 2-component GMM.
    results = comparator._gmm_test_split(data, samples, conditions, 'cpu')
    assert len(results['GMM_chi2_pvalue']) == 1
    assert np.isnan(results['GMM_chi2_pvalue'][0])


def test_gmm_test_split_two_components():
    config = Config(BASIC_CONFIG)
    comparator = TranscriptComparator(config, MockWorker())

    rand_gen = np.random.default_rng(seed=42)
    data = torch.tensor(rand_gen.standard_normal((1, 100, 2)))
    data[:, 50:, :] += 3

    samples = torch.tensor([samp
                            for samp in [0, 1, 2, 3]
                            for _ in range(25)])
    conditions = torch.tensor([cond
                               for cond in [0, 1]
                               for _ in range(50)])

    results = comparator._gmm_test_split(data, samples, conditions, 'cpu')
    assert len(results['GMM_chi2_pvalue']) == 1
    assert results['GMM_chi2_pvalue'] < 0.01


def get_float(value):
    return round(float(value), 3)

