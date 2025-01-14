import torch
import numpy as np

from nanocompore.comparisons import TranscriptComparator
from nanocompore.config import Config
from nanocompore.kmer import KmerData

from tests.common import BASIC_CONFIG


def test_add_shift_stats():
    config = Config(BASIC_CONFIG)
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
    conditions = torch.tensor([[0, 1, 0, 1],
                               [1, 0, 1, 0]])

    comparator = TranscriptComparator(config)

    results = {}
    comparator._add_shift_stats(results,
                                data,
                                conditions,
                                'cpu')

    assert get_float(results['c1_mean_intensity'][0]) == 1.3
    assert get_float(results['c1_mean_intensity'][1]) == 0.9
    assert get_float(results['c1_mean_dwell'][0]) == 0.7
    assert get_float(results['c1_mean_dwell'][1]) == 0.25
    assert get_float(results['c1_mean_motor_dwell'][0]) == 0.5
    assert get_float(results['c1_mean_motor_dwell'][1]) == 0.6

    assert get_float(results['c2_mean_intensity'][0]) == 1.1
    assert get_float(results['c2_mean_intensity'][1]) == 1.2
    assert get_float(results['c2_mean_dwell'][0]) == 0.15
    assert get_float(results['c2_mean_dwell'][1]) == 0.65
    assert get_float(results['c2_mean_motor_dwell'][0]) == 0.65
    assert get_float(results['c2_mean_motor_dwell'][1]) == 0.5

    assert get_float(results['c1_std_intensity'][1]) == 0.3


def test_kmers_to_tensor():
    config = Config(BASIC_CONFIG)
    comparator = TranscriptComparator(config)

    kmers = [
        KmerData('transcript1',
                 3,
                 'AGCAC',
                 np.array(['wt1', 'kd1', 'wt2', 'kd2']),
                 np.array([1, 2, 3, 4]),
                 np.array([97.5, 98.8, 84.3, 113.2]),
                 np.array([3.1, 4.2, 2.1, 1.9]),
                 np.array([0.08, 0.012, 0.4, 0.06]),
                 None,
                 config),
        # we intentionally skip a position to check
        # if the function handles gaps
        KmerData('transcript1',
                 5,
                 'AGCAC',
                 np.array(['wt1', 'kd1', 'wt2', 'kd2']),
                 np.array([1, 2, 3, 4]),
                 np.array([98.5, 99.8, 85.3, 114.2]),
                 np.array([4.1, 5.2, 3.1, 2.9]),
                 np.array([0.09, 0.013, 0.5, 0.07]),
                 None,
                 config),
        KmerData('transcript1',
                 14,
                 'AGCAC',
                 np.array(['wt2', 'kd2', 'wt1', 'kd1']),
                 # note that the reads are in different order here
                 np.array([5, 6, 1, 2]),
                 np.array([96.5, 97.8, 82.3, 111.2]),
                 np.array([2.1, 3.2, 1.1, 0.9]),
                 np.array([0.07, 0.01, 0.3, 0.03]),
                 None,
                 config)
    ]

    data, samples, conditions, positions = comparator._kmers_to_tensor(kmers, 'cpu')

    assert data.shape == (3, 6, 3)

    expected = torch.tensor([
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
    # Apparently torch doesn't have a function to compare
    # two tensors that handles nans...
    assert torch.equal(
        torch.where(data.isnan(), 0, data),
        torch.where(expected.isnan(), 0, expected))


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


def get_float(value):
    return round(float(value), 3)
