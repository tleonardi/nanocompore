import os

from pathlib import Path

import numpy as np

from nanocompore.config import Config
from nanocompore.kmer import KmerData
from nanocompore.preprocessing import Preprocessor
from nanocompore.preprocessing import Uncalled4Preprocessor
from tests.common import BASIC_CONFIG, cwd


def test_get_references_from_bams():
    config = Config(BASIC_CONFIG)
    try:
        preprocessor = Preprocessor(config)
        references = preprocessor._get_references_from_bams()
        assert references == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                          'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}
    finally:
        Path(BASIC_CONFIG['preprocessing_db']).unlink()


def test_get_reads_invalid_kmer_ratio_uncalled4():
    config = Config(BASIC_CONFIG)
    try:
        preprocessor = Uncalled4Preprocessor(config)
        kmer_1 = KmerData('transcript1',
                          1,
                          'AGCAC',
                          np.array(['WT', 'KD', 'WT', 'KD']),
                          np.array([1, 2, 3, 4]),
                          np.array([97.5, 98.8, 84.3, 113.2]),
                          np.array([3.1, 4.2, 2.1, 1.9]),
                          np.array([0.08, 0.012, 0.4, 0.06]),
                          # Validity information is not available
                          # from Uncalled4 and Remora
                          None,
                          config)
        kmer_10 = KmerData('transcript1',
                           10,
                           'AGCAC',
                           np.array(['WT', 'KD', 'WT', 'KD']),
                           np.array([1, 2, 3, 4]),
                           np.array([97.5, 98.8, 84.3, 113.2]),
                           np.array([3.1, 4.2, 2.1, 1.9]),
                           np.array([0.08, 0.012, 0.4, 0.06]),
                          # Validity information is not available
                          # from Uncalled4 and Remora
                           None,
                           config)

        invalid_ratios = preprocessor._get_reads_invalid_kmer_ratio([kmer_1, kmer_10])
        assert invalid_ratios == {1: 0.8,
                                  2: 0.8,
                                  3: 0.8,
                                  4: 0.8}
    finally:
        Path(BASIC_CONFIG['preprocessing_db']).unlink()

