import copy
import sqlite3

from pathlib import Path

import numpy as np

from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import REMORA_MEASUREMENT_TYPE
from nanocompore.common import UNCALLED4_MEASUREMENT_TYPE
from nanocompore.config import Config
from nanocompore.kmer import KmerData
from nanocompore.preprocessing import RemoraPreprocessor
from tests.common import BASIC_CONFIG


# def test_get_reads_invalid_kmer_ratio_uncalled4():
#     yaml = copy.deepcopy(BASIC_CONFIG)
#     yaml['preprocessing_db'] += '.tmp'
#     config = Config(yaml)
#     try:
#         preprocessor = Preprocessor(config)
#         kmer_1 = KmerData('transcript1',
#                           1,
#                           'AGCAC',
#                           np.array(['WT', 'KD', 'WT', 'KD']),
#                           np.array([1, 2, 3, 4]),
#                           np.array([97.5, 98.8, 84.3, 113.2]),
#                           np.array([3.1, 4.2, 2.1, 1.9]),
#                           np.array([0.08, 0.012, 0.4, 0.06]),
#                           # Validity information is not available
#                           # from Uncalled4 and Remora
#                           None,
#                           config)
#         kmer_10 = KmerData('transcript1',
#                            10,
#                            'AGCAC',
#                            np.array(['WT', 'KD', 'WT', 'KD']),
#                            np.array([1, 2, 3, 4]),
#                            np.array([97.5, 98.8, 84.3, 113.2]),
#                            np.array([3.1, 4.2, 2.1, 1.9]),
#                            np.array([0.08, 0.012, 0.4, 0.06]),
#                           # Validity information is not available
#                           # from Uncalled4 and Remora
#                            None,
#                            config)

#         invalid_ratios = preprocessor._get_reads_invalid_kmer_ratio([kmer_1, kmer_10])
#         assert invalid_ratios == {1: 0.8,
#                                   2: 0.8,
#                                   3: 0.8,
#                                   4: 0.8}
#     finally:
#         Path(config.get_preprocessing_db()).unlink()


def test_remora_preprocessing():
    """
    Integration test for RemoraPreprocessor
    """
    yaml = copy.deepcopy(BASIC_CONFIG)
    yaml['preprocessing_db'] += '.tmp'
    yaml['resquiggler'] = 'remora'
    config = Config(yaml)

    preprocessor = RemoraPreprocessor(config)

    try:
        preprocessor()

        conn = sqlite3.connect(config.get_preprocessing_db())

        num_transcripts = conn.execute('SELECT COUNT(*) FROM transcripts').fetchone()[0]
        assert num_transcripts == 2

        kmer_data = conn.execute('SELECT * FROM kmer_data').fetchall()
        assert len(kmer_data) > 0

        # Validate that we have data for two transcripts
        assert {1, 2} == {row[0] for row in kmer_data}

        # Validate that all reads are included in the
        # preprocessing database.
        read_ids = {int(read_id)
                    for row in kmer_data
                    for read_id in np.frombuffer(row[4], dtype=READ_ID_TYPE)}
        query = f'''SELECT read
                    FROM reads
                    WHERE id IN ({",".join(["?" for _ in range(len(read_ids))])})'''
        read_uuids = conn.execute(query, tuple(read_ids)).fetchall()
        read_uuids = {row[0] for row in read_uuids}

        assert read_uuids == {
                # kd1
                '3f46f499-8ce4-4817-8177-8ad61b784f27',
                '73d62df4-f04a-4207-a4bc-7b9739b3c3b2',
                'b7bc9a36-318e-4be2-a90f-74a5aa6439bf',
                '65db1ced-cbf5-40e4-880b-9f202c715804',
                '6384da30-e6cc-4421-b93a-3aca69d09857',
                '6c80b542-3176-49c8-b6df-a204499cc0dd',
                # kd2
                'ac486e16-15be-47a8-902c-2cfa2887c534',
                '797fd991-570e-42d4-8292-0a7557b192d7',
                '4e1ad358-ec2b-40b4-8e9a-54db28a40551',
                '5490c55c-906d-4137-aaf8-1bec590443db',
                'c7e58c3b-5007-4e6f-aa68-eee02df989ab',
                '140112e7-222a-45e8-bbd2-37d1a724fb09',
                # wt1
                'a4395b0d-dd3b-48e3-8afb-4085374b1147',
                'f9733448-6e6b-47ba-9501-01eda2f5ea26',
                '6f5e3b2e-f27b-47ef-b3c6-2ab4fdefd20a',
                'bfaf670f-2c38-4984-93fc-5f704ce7cb5d',
                '07cf4a88-6083-462b-a010-1ef095a4e457',
                'b9f46baf-d038-46be-b55a-7f3a82716d59',
                # wt2
                '2da07406-70c2-40a1-835a-6a7a2c914d49',
                '54fc1d38-5e3d-4d77-a717-2d41b4785af6',
                '3cfa90d1-7dfb-4398-a224-c75a3ab99873',
                '598b375f-cdfe-41c8-aca2-78ffc0adbaf7',
                'a75173bb-74ca-4b38-a579-aeddee90518e',
                '8a3709ea-cd84-4b19-a33e-83198f9f323f'}


        get_actb_id = 'SELECT id FROM transcripts WHERE name LIKE "%ENST00000674681.1%"'
        transcript_id = conn.execute(get_actb_id).fetchone()[0]

        kmer_301 = [row
                    for row in kmer_data
                    if row[1] == 301 and row[0] == transcript_id][0]

        # Validate that the kmer sequence is correct
        assert kmer_301[2] == 'A'

        intensities = np.frombuffer(kmer_301[5], dtype=REMORA_MEASUREMENT_TYPE)
        assert round(float(np.mean(intensities)), 3) == 0.8
    finally:
        Path(config.get_preprocessing_db()).unlink()

