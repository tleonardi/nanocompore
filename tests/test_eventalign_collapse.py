import os
import sqlite3

from pathlib import Path

import numpy as np

from nanocompore.eventalign_collapse import EventalignCollapser
from nanocompore.common import MEASUREMENTS_TYPE

from tests.common import cwd


def test_eventalign_collapse():
    """
    An integration tests that runs the whole
    eventalign_collapse process.
    """
    tsv_in = os.path.join(cwd, 'tests/fixtures/kd1_eventalign.tsv')
    ref = os.path.join(cwd, 'tests/fixtures/test_reference.fa')
    db_out = 'collapsed.sqlite'

    collapser = EventalignCollapser(tsv_in, ref, db_out, 2, '.')

    try:
        collapser()

        conn = sqlite3.connect(db_out)
        num_transcripts = conn.execute('SELECT COUNT(*) FROM transcripts').fetchone()[0]
        assert num_transcripts == 2

        signal_data = conn.execute('SELECT * FROM signal_data').fetchall()

        # Validate that we have data for two transcripts
        assert {1, 2} == {row[0] for row in signal_data}

        # Validate that we have four reads
        read_ids = [row[1] for row in signal_data]
        assert len(read_ids) == 4

        query = f'''SELECT read
                    FROM reads
                    WHERE id IN ({",".join(["?" for _ in range(len(read_ids))])})'''
        read_uuids = conn.execute(query, read_ids).fetchall()
        read_uuids = {row[0] for row in read_uuids}

        assert read_uuids == {'3f46f499-8ce4-4817-8177-8ad61b784f27',
                              '73d62df4-f04a-4207-a4bc-7b9739b3c3b2',
                              'b7bc9a36-318e-4be2-a90f-74a5aa6439bf',
                              '65db1ced-cbf5-40e4-880b-9f202c715804'}

        # Validate the measurements for the first two
        # positions of the first read.
        intensities = np.frombuffer(signal_data[0][2], dtype=MEASUREMENTS_TYPE)
        assert round(float(intensities[301]), 3) == 117.218
        assert round(float(intensities[302]), 3) == 71.906

        dwells = np.frombuffer(signal_data[0][3], dtype=MEASUREMENTS_TYPE)
        assert round(float(dwells[301]), 3) == 0.008
        assert round(float(dwells[302]), 3) == 0.003

        # Validate the calculation of read invalid
        # kmer ratio.
        reads_ratios = {read: invalid_ratio
                        for read, _, invalid_ratio in conn.execute('SELECT * FROM reads').fetchall()}
        # Read 3f46f499-8ce4-4817-8177-8ad61b784f27
        # is for transcript:
        # ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|
        # We have three positions in the eventalign tsv
        # with positions 302 and 303 being valid.
        # Position 301 has two events, one of which is
        # NNNNN, so the position is invalid.
        # The expected invalid ratio should be
        # the (alignment_length - num_valid)/alignment_length
        assert reads_ratios['3f46f499-8ce4-4817-8177-8ad61b784f27'] == (3 - 2)/3
    finally:
        # Make sure to delete the sqlite database
        Path(db_out).unlink()


def test_read_data():
    tsv_in = os.path.join(cwd, 'tests/fixtures/kd1_eventalign.tsv')
    ref = os.path.join(cwd, 'tests/fixtures/test_reference.fa')
    db_out = 'collapsed.sqlite'

    collapser = EventalignCollapser(tsv_in, ref, db_out, 2, '.')

    fstart = 161
    fend = 6649

    with open(tsv_in) as f:
        data = collapser._read_data(f, fstart, fend)

    read, kmers = data[0]

    assert read == '3f46f499-8ce4-4817-8177-8ad61b784f27'
    assert len(kmers) == 3
    assert kmers[0].pos == 301

    # Test that the two events for the position
    # are collapsed:
    assert kmers[0].dwell == 0.00498 + 0.00332
    # The median intensity should be the median
    # of all sample values from the two events.
    assert kmers[0].median_intensity == 117.218
    # The first event is valid, but the second
    # is not, so the collapsed kmer should be
    # invalid as well.
    assert not kmers[0].valid


def test_process_ref():
    tsv_in = os.path.join(cwd, 'tests/fixtures/kd1_eventalign.tsv')
    ref = os.path.join(cwd, 'tests/fixtures/test_reference.fa')
    db_out = 'collapsed.sqlite'

    collapser = EventalignCollapser(tsv_in, ref, db_out, 2, '.')

    fstart = 161
    fend = 6649

    with open(tsv_in) as f:
        data = collapser._read_data(f, fstart, fend)

    ref_id = 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|'
    ref_len = 2554

    intensity, dwell, read_invalid_ratios = collapser._process_ref(ref_id, data, ref_len)

    assert intensity.shape == (3, 2554)
    assert round(float(intensity[0, 301]), 3) == 117.218

    assert dwell.shape == (3, 2554)
    assert round(float(dwell[0, 301]), 4) == 0.00498 + 0.00332

    assert len(read_invalid_ratios) == 3
    assert round(read_invalid_ratios['3f46f499-8ce4-4817-8177-8ad61b784f27'], 3) == 0.333

