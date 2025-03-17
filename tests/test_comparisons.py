import copy

import torch
import numpy as np
import pandas as pd

from nanocompore.comparisons import TranscriptComparator
from nanocompore.comparisons import calculate_lor
from nanocompore.comparisons import crosstab
from nanocompore.config import Config

from tests.common import BASIC_CONFIG


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
                               

    comparator = TranscriptComparator(config)

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

    assert get_float(results['c1_std_intensity'][1]) == 0.1


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

    comparator = TranscriptComparator(config)

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

    assert get_float(results['c1_std_intensity'][1]) == 0.1


def test_nonparametric_test_KS():
    config = Config(BASIC_CONFIG)
    comparator = TranscriptComparator(config)
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
    comparator = TranscriptComparator(config)
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
    contingency = np.array([[13, 11],
                            [12, 10]])
    assert calculate_lor(contingency) == -0.015


def test_crosstab():
    table = crosstab(np.array([0, 1, 0, 1, 0, 1]),
                     np.array([0, 0, 0, 1, 1, 1]))
    assert table[0, 0] == 2
    assert table[0, 1] == 1
    assert table[1, 0] == 1
    assert table[1, 1] == 2


def test_get_cluster_counts():
    config = Config(BASIC_CONFIG)
    # Condition 0 is the knockdown (depleted condition).
    # Since most of the points for the non-depleted condition (1)
    # are found in cluster 1, then we can guess that the
    # modification cluster is 1.
    conditions = torch.tensor([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
    samples = torch.tensor([0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3])
    predictions = torch.tensor([0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1])

    contingency = crosstab(conditions, predictions)

    comparator = TranscriptComparator(config)


    cluster_counts = comparator._get_cluster_counts(contingency, samples, conditions, predictions)

    print(cluster_counts)
    assert cluster_counts == {'kd1_mod': 0, 'kd1_unmod': 2,
                              'kd2_mod': 1, 'kd2_unmod': 2,
                              'wt1_mod': 1, 'wt1_unmod': 1,
                              'wt2_mod': 3, 'wt2_unmod': 1}


def get_float(value):
    return round(float(value), 3)

