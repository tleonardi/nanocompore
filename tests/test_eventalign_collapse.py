import os
import sqlite3

from pathlib import Path

import numpy as np

from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.common import Kit
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.kmer import KmerData

from tests.common import BASIC_CONFIG, cwd


def test_eventalign_collapse():
    """
    An integration tests that runs the whole
    eventalign_collapse process.
    """
    tsv_in = os.path.join(cwd, 'tests/fixtures/kd1_eventalign.tsv')
    ref = os.path.join(cwd, 'tests/fixtures/test_reference.fa')
    db_out = 'collapsed.sqlite'

    collapser = EventalignCollapser(tsv_in, ref, db_out, Kit.RNA002, 2)

    try:
        collapser()

        conn = sqlite3.connect(db_out)
        num_transcripts = conn.execute('SELECT COUNT(*) FROM transcripts').fetchone()[0]
        assert num_transcripts == 2

        kmer_data = conn.execute('SELECT * FROM kmer_data').fetchall()
        # Validate that we have data for two transcripts
        assert {1, 2} == {row[0] for row in kmer_data}

        # Validate that we have data for three positions and the
        # center of kmers is correct. The positions in the
        # eventalign tsv are 301, 302 and 303, but eventalign
        # uses the position of the first base in the k-mer.
        # With RNA002, we have 5mers and the center (most-influential)
        # position for that chemistry kit is the 4rd one.
        assert {304, 305, 306} == {row[1] for row in kmer_data}

        kmer_304 = [row
                    for row in kmer_data
                    if row[1] == 304 and row[0] == 1][0]

        # Validate that the kmer sequence is correct
        assert kmer_304[2] == 'AGCAC'

        read_ids = np.frombuffer(kmer_304[3], dtype=READ_ID_TYPE)
        assert len(read_ids) == 3

        query = f'''SELECT read
                    FROM reads
                    WHERE id IN ({",".join(["?" for _ in range(3)])})'''
        read_uuids = conn.execute(query, read_ids.tolist()).fetchall()
        read_uuids = {row[0] for row in read_uuids}

        assert read_uuids == {'3f46f499-8ce4-4817-8177-8ad61b784f27',
                              '73d62df4-f04a-4207-a4bc-7b9739b3c3b2',
                              'b7bc9a36-318e-4be2-a90f-74a5aa6439bf'}

        # Validate the median intensity for position
        # 304 on the first read (which includes two events).
        intensities = np.frombuffer(kmer_304[4], dtype=EVENTALIGN_MEASUREMENT_TYPE)[0]
        assert round(float(np.median(intensities)), 3) == 117.218

        # Validate the median intensity for position
        # 305 on the first read (which includes one event).
        kmer_305 = [row
                    for row in kmer_data
                    if row[1] == 305 and row[0] == 1][0]
        intensities = np.frombuffer(kmer_305[4], dtype=EVENTALIGN_MEASUREMENT_TYPE)[0]
        assert round(float(np.median(intensities)), 3) == 71.906

        # Validate the MAD calcuation
        intensity_mads = np.frombuffer(kmer_304[5], dtype=EVENTALIGN_MEASUREMENT_TYPE)
        assert round(float(intensity_mads[0]), 3) == 4.723

        dwells = np.frombuffer(kmer_304[6], dtype=EVENTALIGN_MEASUREMENT_TYPE)
        assert round(float(dwells[0]), 3) == 0.008

    finally:
        # Make sure to delete the sqlite database
        Path(db_out).unlink()


def test_get_reads_invalid_kmer_ratio():
    collapser = EventalignCollapser(None,
                                    None,
                                    None,
                                    Kit.RNA002,
                                    2)
    kmer_1 = KmerData('transcript1',
                      1,
                      'AGCAC',
                      np.array(['WT', 'KD', 'WT', 'KD']),
                      np.array([1, 2, 3, 4]),
                      np.array([97.5, 98.8, 84.3, 113.2]),
                      np.array([3.1, 4.2, 2.1, 1.9]),
                      np.array([0.08, 0.012, 0.4, 0.06]),
                      np.array([True, False, True, False]),
                      None)
    kmer_10 = KmerData('transcript1',
                       10,
                       'AGCAC',
                       np.array(['WT', 'KD', 'WT', 'KD']),
                       np.array([1, 2, 3, 4]),
                       np.array([97.5, 98.8, 84.3, 113.2]),
                       np.array([3.1, 4.2, 2.1, 1.9]),
                       np.array([0.08, 0.012, 0.4, 0.06]),
                       np.array([True, False, False, True]),
                       None)
    invalid_ratios = collapser._get_reads_invalid_kmer_ratio([kmer_1, kmer_10], 10)
    assert invalid_ratios == {1: (10 - 2)/10,
                              2: 1,
                              3: (10 - 1)/10,
                              4: (10 - 1)/10}

