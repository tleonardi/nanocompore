import sqlite3
import json
from contextlib import closing

import numpy as np
import pandas as pd

from nanocompore.database import DB_METADATA_READ_ID_TYPE_KEY
from nanocompore.database import DB_METADATA_SAMPLE_ID_TYPE_KEY
from nanocompore.database import DB_METADATA_MEASUREMENT_TYPE_KEY
from nanocompore.database import DB_METADATA_SAMPLE_LABELS_KEY
from nanocompore.database import DB_METADATA_CONDITION_SAMPLES_KEY

SAMPLE_TO_CONDITION_KEY = 'sample_to_condition'


GET_KMER_SQL = """
SELECT *
FROM kmer_data
WHERE transcript_id = ?
  AND pos = ?
"""

GET_READS_SQL = """
SELECT pos, samples, reads, {}
FROM kmer_data
WHERE transcript_id = ?
"""

GET_MINMAX_POS = """
SELECT MIN(pos), MAX(pos)
FROM kmer_data
WHERE transcript_id = ?
"""

GET_TRANSCRIPT_ID = """
SELECT id FROM transcripts WHERE reference = ?
"""


def get_pos(db, reference_id, pos):
    """
    Returns the signal data for a specific position
    of the given reference transcript from all reads.

    Returns: (signal_dataframe, kmer)
             signal_dataframe contains the columns:
                - condition: the read's condition label
                - sample:    the read's sample label
                - read:      id of the read
                - intensity: current intensity
                - std:       std of the intensity
                - dwell:     dwell time for the kmer
    """
    metadata = get_metadata(db)

    with closing(sqlite3.connect(db)) as conn,\
         closing(conn.cursor()) as cursor:
        transcript_id = cursor.execute(GET_TRANSCRIPT_ID,
                                       (reference_id,)).fetchone()[0]
        res = cursor.execute(GET_KMER_SQL, (transcript_id, pos))
        row = res.fetchone()

    kmer = row[2]
    sample = _parse_samples(row[3], metadata)
    return pd.DataFrame({
        'condition': [metadata[SAMPLE_TO_CONDITION_KEY][samp] for samp in sample],
        'sample': sample,
        'read': _parse_reads(row[4], metadata),
        'intensity': _parse_measurements(row[5], metadata),
        'std': _parse_measurements(row[6], metadata),
        'dwell': _parse_measurements(row[7], metadata),
    }), kmer


def get_read(db,
             reference_id,
             read_id,
             positions=None,
             variables=('dwell', 'intensity')):
    """
    Get a numpy array with the signal data for a specific read from a kmer data DB.

    By default it returns intensity and dwell for all positions.
    Passing a collection of positions filtering the results.
    Possible variables are ("intensity", "intensity_std", "dwell").

    Returns: (signal_data, sample, condition)

             signal_data is a 2D array with shape (positions, variables).
             sample is the sample label for the read
             conditions is the condition label for the read

    Note: the positions would always span from 0 to the maximum position
    with data available. If no signal data is available for a position
    it will contain a np.nan value.
    """
    signal_data, _, samples =  get_reads(db,
                                         reference_id,
                                         reads=[read_id],
                                         positions=positions,
                                         variables=variables)
    return signal_data[0], samples[0]


def get_reads(db,
              reference_id,
              positions=None,
              reads=None,
              variables=('dwell', 'intensity')):
    """
    Get a numpy array with read data for specific reference_id from
    a kmer_data sqlite DB.

    By default it returns intensity and dwell for all positions and reads.
    Passing a collection of positions and reads can be used for filtering
    the results.
    Possible variables are ("intensity", "intensity_std", "dwell").


    Returns: (signal_data, reads, samples, conditions)

             signal_data is a 3D array with shape (reads, positions, variables).
             reads is an array of read uuids
             samples is an array of sample labels
             conditions is an array of condition labels

    Note: the positions would always span from 0 to the maximum position
    with data available. If no signal data is available for a position
    on a read, the array will contain np.nan at the corresponding indices.
    """
    if positions:
        positions = set(positions)
    if reads:
        reads = set(reads)

    for var in variables:
        if var not in ['intensity', 'intensity_std', 'dwell']:
            raise ValueError(f"Variable {var} not supported. You can use 'intensity', 'intensity_std' and 'dwell'.")

    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    metadata = get_metadata(db)

    transcript_id = cursor.execute(GET_TRANSCRIPT_ID,
                                   (reference_id,)).fetchone()[0]

    if positions:
        ref_len = max(positions)
    else:
        _, ref_len = cursor.execute(GET_MINMAX_POS, (transcript_id,)).fetchone()

    res = cursor.execute(GET_READS_SQL.format(', '.join(variables)),
                         (transcript_id,))

    data = [_map_db_kmer_row(row, variables, metadata)
            for row in res.fetchall()
            if positions is None or row[0] in positions]

    # Get all reads for which there's data for this transcript
    all_read_ids = list({read
                         for row in data
                         for read in row['reads']})

    # Map the integer ids of the reads to their uuids
    ids_to_reads = _get_read_uuids(cursor, all_read_ids)

    # If the user requests specific reads,
    # we filter all reads to include only them.
    if reads:
        all_read_ids = [read_id
                        for read_id in all_read_ids
                        if ids_to_reads[read_id] in reads]
    nreads = len(all_read_ids)

    # Fix indices in the tensor for each read and variable
    read_indices = dict(zip(all_read_ids, range(nreads)))
    var_indices = dict(zip(variables, range(len(variables))))

    tensor = np.empty((nreads, ref_len, len(variables)))
    tensor.fill(np.nan)
    res_reads = np.empty(nreads, dtype=object)
    res_samples = np.empty(nreads, dtype=object)
    res_conditions = np.empty(nreads, dtype=object)

    for kmer in data:
        pos = kmer['pos']
        kmer_samples = kmer['samples']
        kmer_conditions = kmer['conditions']
        kmer_read_ids = kmer['reads']
        for i, read_id in enumerate(kmer_read_ids):
            read = ids_to_reads[read_id]
            if reads and read not in reads:
                continue
            read_indx = read_indices[read_id]
            res_reads[read_indx] = read
            res_samples[read_indx] = kmer_samples[i]
            res_conditions[read_indx] = kmer_conditions[i]
            for var in variables:
                var_indx = var_indices[var]
                tensor[read_indx, pos-1, var_indx] = kmer[var][i]

    return tensor, res_reads, res_samples, res_conditions


def get_metadata(db):
    with closing(sqlite3.connect(db)) as conn,\
         closing(conn.cursor()) as cursor:
        query = "SELECT key, value FROM metadata"
        metadata = {k: v
                    for k, v in cursor.execute(query).fetchall()}
        metadata[DB_METADATA_CONDITION_SAMPLES_KEY] = json.loads(metadata[DB_METADATA_CONDITION_SAMPLES_KEY])
        metadata[SAMPLE_TO_CONDITION_KEY] = {s: c
                                             for c, samples in metadata[DB_METADATA_CONDITION_SAMPLES_KEY].items()
                                             for s in samples}
        return metadata


def _map_db_kmer_row(row, variables, metadata):
    pos = row[0]
    samples = _parse_samples(row[1], metadata)
    read_ids = _parse_reads(row[2], metadata)
    result = {
        'pos': pos,
        'samples': samples,
        'conditions': [metadata[SAMPLE_TO_CONDITION_KEY][samp] for samp in samples],
        'reads': read_ids,
    }
    for i, var in enumerate(variables):
        result[var] = _parse_measurements(row[3+i], metadata)
    return result


def _parse_reads(binary, metadata):
    read_id_type = getattr(np, metadata[DB_METADATA_READ_ID_TYPE_KEY])
    return np.frombuffer(binary, dtype=read_id_type)


def _parse_samples(binary, metadata):
    sample_id_type = getattr(np, metadata[DB_METADATA_SAMPLE_ID_TYPE_KEY])
    sample_labels = metadata[DB_METADATA_SAMPLE_LABELS_KEY].split(',')
    return np.array([sample_labels[i]
                    for i in np.frombuffer(binary, dtype=sample_id_type)])


def _parse_measurements(binary, metadata):
    measurement_type = getattr(np, metadata[DB_METADATA_MEASUREMENT_TYPE_KEY])
    return np.frombuffer(binary, dtype=measurement_type)


def _get_read_uuids(cursor, reads_ids):
    query = "SELECT id, read FROM reads WHERE id IN ({})"
    query = query.format(','.join(map(str, reads_ids)))
    mappings = {}
    for ids in _chunk(reads_ids, 1000):
        for row in cursor.execute(query).fetchall():
            mappings[row[0]] = row[1]

    return mappings


def _chunk(li, n):
    for i in range(0, len(li), n):
        yield li[i:i + n]

