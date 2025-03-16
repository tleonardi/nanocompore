import sqlite3
from typing import Optional
from contextlib import closing

import numpy as np
import pandas as pd
from jaxtyping import Float

from nanocompore.common import MEASUREMENTS_TYPE
from nanocompore.database import PreprocessingDB


def get_references(db: str, has_data=True) -> list[str]:
    """
    Returns a list of all references in the given
    database.

    Parameters
    ----------
    db : str
        Path to the SQLite database produced by the preprocessing
        command of Nanocompore.
    has_data : bool
        If True (default) will return only
        references for which there are data.

    Returns
    -------
    list
        List of transcript reference id strings.
    """
    if has_data:
        transcripts = PreprocessingDB(db).get_references_with_data()
    else:
        transcripts = PreprocessingDB(db).get_references()
    return [t.ref_id for t in transcripts]


def get_reads(db: str,
              reference_id: str,
              selected_reads: Optional[list[str]]=None,
              variables: Optional[tuple[str, ...]]=('dwell', 'intensity')
    ) -> tuple[
            Float[np.ndarray, "reads positions variables"],
            list[str]]:
    """
    Get a numpy array with read data for specific reference_id from
    a preprocessing database produced by eventalign_collapse or
    remora_resquiggle.

    Parameters
    ----------
    db : str
        Path to the SQLite database file.
    reference_id : str
        ID for a reference sequence (transcript).
    selected_reads : Optional[list[str]]
        Optional list of UUIDs of the reads for which to get data.

    Returns
    -------
    tuple[Float[np.ndarray, ["reads positions variables"]],
          list[str]]

        A tuple with (signal_data, read_ids)
        - signal_data is a 3D array with shape (reads, positions, variables).
          In the variables dimension 0=intensity, 1=dwell time.
        - reads is a list of read ids (qname).
    """
    with closing(sqlite3.connect(db)) as conn,\
         closing(conn.cursor()) as cursor:
        rows = cursor.execute("""
            SELECT intensity, dwell, read
            FROM signal_data sd
            JOIN transcripts t ON t.id = sd.transcript_id
            JOIN reads r ON r.id = sd.read_id
            WHERE t.name = ?""", (reference_id,)).fetchall()
        intensity = np.array([np.frombuffer(row[0], dtype=MEASUREMENTS_TYPE)
                              for row in rows])
        dwell = np.array([np.frombuffer(row[1], dtype=MEASUREMENTS_TYPE)
                          for row in rows])
        reads = [row[2] for row in rows]
        signal_data = np.stack([intensity, dwell], axis=2)
        if selected_reads:
            read_set = set(selected_reads)
            mask = np.array([read in selected_reads for read in reads])
            filtered_reads = [read for read in reads if read in selected_reads]
            return signal_data[mask], filtered_reads
        return signal_data, reads


def get_read(db: str, read_id: str) -> Float[np.ndarray, "reads positions variables"]:
    """
    Get a numpy array with the signal data for a specific read from
    a preprocessing database produced by eventalign_collapse or
    remora_resquiggle.

    Parameters
    ----------
    db : str
        Path to the SQLite database file.
    read_id : str
        ID of the read for which to get data.

    Returns
    -------
    Float[np.ndarray, "positions variables"]
        Sinal data is a 2D numpy array with shape
        (positions, variables). Variables
        will be 0=intensity and 1=dwell time.
    """
    with closing(sqlite3.connect(db)) as conn,\
         closing(conn.cursor()) as cursor:
        row = cursor.execute("""
            SELECT intensity, dwell
            FROM signal_data sd
            JOIN reads r ON r.id = sd.read_id
            WHERE r.read = ?
            LIMIT 1""", (read_id,)).fetchone()
        intensity = np.frombuffer(row[0], dtype=MEASUREMENTS_TYPE)
        dwell = np.frombuffer(row[1], dtype=MEASUREMENTS_TYPE)
        return np.stack([intensity, dwell], axis=1)


def get_pos(db: str, reference_id: str, pos: int) -> pd.DataFrame:
    """
    Get data for a position.

    Returns the signal data for a specific position
    of the given reference transcript from all reads.

    Parameters
    ----------
    db : str
        Path to the SQLite database produced by the preprocessing
        command of Nanocompore.
    reference_id : str
        ID for a reference sequence (transcript).
    pos : int
        Position on the transcript for which to get data.

    Returns
    -------
    pandas.DataFrame
        Where the DataFrame contains the following columns:
            - read:      id of the read
            - intensity: current intensity
            - dwell:     dwell time for the kmer
    """
    reads_data, reads = get_reads(db, reference_id)
    pos_data = reads_data[:, pos, :]
    valid = ~np.isnan(pos_data[:, 0])
    return pd.DataFrame({'read': np.array(reads)[valid],
                         'intensity': pos_data[valid, 0],
                         'dwell': pos_data[valid, 1]})


def get_metadata(db: str) -> dict[str, str]:
    """
    Returns the metadata from the given database.

    The metadata contains information such as input files,
    resquiggler used, and data types for the binary encoded fields.

    Parameters
    ----------
    db : str
        Path to the SQLite database produced by the preprocessing
        command of Nanocompore.

    Returns
    -------
    dict
        Dictionary containing the metadata
    """
    with closing(sqlite3.connect(db)) as conn,\
         closing(conn.cursor()) as cursor:
        query = "SELECT key, value FROM metadata"
        return {k: v for k, v in cursor.execute(query).fetchall()}

