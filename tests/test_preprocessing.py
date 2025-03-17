import copy
import sqlite3

from pathlib import Path

import numpy as np

from nanocompore.common import READ_ID_TYPE
from nanocompore.common import MEASUREMENTS_TYPE
from nanocompore.common import Kit
from nanocompore.config import Config
from nanocompore.preprocessing import RemoraPreprocessor
from tests.common import BASIC_CONFIG


def test_remora_preprocessing():
    """
    Integration test for RemoraPreprocessor
    """
    output = 'kd1_remora.sqlite'
    preprocessor = RemoraPreprocessor(BASIC_CONFIG['fasta'],
                                      'tests/fixtures/kd1.pod5',
                                      'tests/fixtures/kd1.bam',
                                      output,
                                      Kit['RNA002'],
                                      5000,
                                      2)

    try:
        preprocessor()

        conn = sqlite3.connect(output)

        num_transcripts = conn.execute('SELECT COUNT(*) FROM transcripts').fetchone()[0]
        assert num_transcripts == 2

        signal_data = conn.execute('SELECT * FROM signal_data').fetchall()
        assert len(signal_data) > 0

        # Validate that we have data for two transcripts
        assert {1, 2} == {row[0] for row in signal_data}

        # Validate that all reads are included in the
        # preprocessing database.
        read_ids = {row[1] for row in signal_data}
        query = f'''SELECT read
                    FROM reads
                    WHERE id IN ({",".join(["?" for _ in range(len(read_ids))])})'''
        read_uuids = conn.execute(query, tuple(read_ids)).fetchall()
        read_uuids = {row[0] for row in read_uuids}

        assert read_uuids == {
                '3f46f499-8ce4-4817-8177-8ad61b784f27',
                '73d62df4-f04a-4207-a4bc-7b9739b3c3b2',
                'b7bc9a36-318e-4be2-a90f-74a5aa6439bf',
                '65db1ced-cbf5-40e4-880b-9f202c715804',
                '6384da30-e6cc-4421-b93a-3aca69d09857',
                '6c80b542-3176-49c8-b6df-a204499cc0dd'}


        get_actb_id = 'SELECT id FROM transcripts WHERE name LIKE "%ENST00000674681.1%"'
        transcript_id = conn.execute(get_actb_id).fetchone()[0]

        # get the intensity for all reads for ACTB
        reads = np.array([np.frombuffer(row[2], dtype=MEASUREMENTS_TYPE)
                          for row in signal_data
                          if row[0] == transcript_id])

        pos301 = reads[:, 301]
        assert round(float(pos301.mean()), 3) == 0.806
    finally:
        Path(output).unlink()

