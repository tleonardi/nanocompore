import sqlite3
import yaml
from typing import Optional, Union
from contextlib import closing

import numpy as np
import pandas as pd
import pysam
from jaxtyping import Float

from nanocompore.common import MEASUREMENTS_TYPE
from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS
from nanocompore.common import UNCALLED4
from nanocompore.common import Kit
from nanocompore.common import get_references_from_bam
from nanocompore.config import Config
from nanocompore.database import PreprocessingDB
from nanocompore.uncalled4 import Uncalled4


def load_config(config_path: str) -> Config:
    """
    Load a configuration file.

    Parameters
    ----------
    config_path : str
        Path to the Nanocompore configuration file.

    Returns
    -------
    Config
        A configuration object.
    """
    with open(config_path, 'rb') as f:
        return Config(yaml.safe_load(f))


def get_references(config: Config, has_data=True) -> list[str]:
    """
    Returns a list of all references found in the
    list of samples defined in the configuration.

    Parameters
    ----------
    config : Config
        Path to a Nanocompore configuration file.
    has_data : bool, default=True
        If True (default) will return only references
        for which there are mapped reads.

    Returns
    -------
    list
        List of transcript reference id strings.

    Examples
    --------
    >>> from nanocompore.api import load_config, get_references
    >>> config = load_config('analysis.yaml')
    >>> get_references(config)
    ['ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|']
    """
    data_files = list(_get_data_files(config).values())
    if config.get_resquiggler() == UNCALLED4:
        return _get_bam_references(data_files, has_data)
    else:
        return _get_db_references(data_files, has_data)


def get_reads(
        config: Config,
        reference_id: str,
        selected_reads: Optional[list[str]]=None
    ) -> tuple[Float[np.ndarray, "reads positions variables"],
               list[str],
               list[str],
               list[str]]:
    """
    Get the data for all reads mapping to the given reference.

    Parameters
    ----------
    config : Config
        Path to a Nanocompore configuration file.
    reference_id : str
        ID for a reference sequence (transcript).
    selected_reads : Optional[list[str]]
        Optional list of UUIDs of the reads for which to get data.
        By default it's set to None and returns all reads.

    Returns
    -------
    tuple[Float[np.ndarray, ["reads positions variables"]],
          list[str],
          list[str],
          list[str]]

        A tuple with (signal_data, reads, samples, conditions)

        - signal_data is a 3D array with shape (reads, positions, variables).
          In the variables dimension 0=intensity, 1=dwell time.
        - reads is a list of read ids (qname).
        - samples is a list of the sample labels (as defined in the config).
        - conditions is a list of the condition labels (as defined in the config).

    Raises
    ------
    KeyError
        If the reference_id is not found in the data sources.

    Examples
    --------
    >>> from nanocompore.api import load_config, get_references
    >>> config = load_config('analysis.yaml')
    >>> get_reads(config, 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')
    """
    data_files = _get_data_files(config)
    if config.get_resquiggler() == UNCALLED4:
        kit = config.get_kit()
        data, reads, samples = _get_bam_reads(data_files.values(),
                                              reference_id,
                                              kit,
                                              selected_reads)
    else:
        data, reads, samples = _get_db_reads(data_files.values(),
                                             reference_id,
                                             selected_reads)
    sample_mapper = np.vectorize(dict(enumerate(data_files)).get)
    samples = sample_mapper(samples)
    condition_mapper = np.vectorize(config.sample_to_condition().get)
    conditions = condition_mapper(samples)
    return data, reads, samples.tolist(), conditions.tolist()


def get_pos(config: Config, reference_id: str, pos: int) -> pd.DataFrame:
    """
    Get the data for a given position for all samples.
    Note that position is a 0-based index of the first
    nucleotide of a k-mer.

    Returns the signal data for a specific position
    of the given reference transcript from all reads.

    Parameters
    ----------
    config : Config
        Path to a Nanocompore configuration file.
    reference_id : str
        ID for a reference sequence (transcript).
    pos : int
        Position on the transcript for which to get data. A 0-based
        index is assumed.

    Returns
    -------
    pandas.DataFrame
        Where the DataFrame contains the following columns:

        - condition  condition label (as defined in the configuration)
        - sample     sample label (as defined in the configuration)
        - read:      id of the read (qname)
        - intensity: current intensity
        - dwell:     dwell time for the kmer

    Examples
    --------
    >>> from nanocompore.api import load_config, get_pos
    >>> config = load_config('analysis.yaml')
    >>> get_pos(config, 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', 532)
       condition sample                                  read  intensity  dwell
    0         WT    WT1  a4395b0d-dd3b-48e3-8afb-4085374b1147     3800.0    7.0
    1         WT    WT1  f9733448-6e6b-47ba-9501-01eda2f5ea26     4865.0  126.0
    2         WT    WT1  6f5e3b2e-f27b-47ef-b3c6-2ab4fdefd20a     3272.0   42.0
    3         WT    WT2  2da07406-70c2-40a1-835a-6a7a2c914d49     6241.0   44.0
    4         WT    WT2  54fc1d38-5e3d-4d77-a717-2d41b4785af6     4047.0    9.0
    5         WT    WT2  3cfa90d1-7dfb-4398-a224-c75a3ab99873     3709.0   70.0
    6         KD    KD1  3f46f499-8ce4-4817-8177-8ad61b784f27     4807.0   57.0
    7         KD    KD1  73d62df4-f04a-4207-a4bc-7b9739b3c3b2     4336.0  132.0
    8         KD    KD1  b7bc9a36-318e-4be2-a90f-74a5aa6439bf     -861.0    7.0
    9         KD    KD2  ac486e16-15be-47a8-902c-2cfa2887c534     2706.0   45.0
    10        KD    KD2  797fd991-570e-42d4-8292-0a7557b192d7     5450.0   24.0
    11        KD    KD2  4e1ad358-ec2b-40b4-8e9a-54db28a40551      206.0   47.0
    """
    data_files = _get_data_files(config)
    sample_mapper = np.vectorize(dict(enumerate(data_files)).get)
    if config.get_resquiggler() == UNCALLED4:
        kit = config.get_kit()
        df = _get_bam_pos(data_files.values(), reference_id, pos, kit)
    else:
        df = _get_db_pos(data_files.values(), reference_id, pos)
    df['sample'] = sample_mapper(df['sample'])
    condition_mapper = np.vectorize(config.sample_to_condition().get)
    df['condition'] = condition_mapper(df['sample'])
    return df.loc[:, ['condition', 'sample', 'read', 'intensity', 'dwell']]


def _get_bam_references(bams: Union[str, list[str]], has_data=True) -> list[str]:
    """
    Returns a list of all references in the given
    BAM.

    Parameters
    ----------
    bams : Union[str, list[str]]
        One on more paths to BAM files produced by Uncalled4.
    has_data : bool, default=True
        If True (default) will return only references
        for which there are mapped reads.

    Returns
    -------
    list
        List of transcript reference id strings.

    Examples
    --------
    >>> from nanocompore.api import get_bam_references
    >>> get_bam_references('wt1.bam')
    ['ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|']
    """
    if isinstance(bams, str):
        return list(get_references_from_bam(bams, has_data).keys())
    return list({ref
                 for bam in bams
                 for ref in get_references_from_bam(bam, has_data).keys()})


def _get_db_references(dbs: Union[str, list[str]], has_data=True) -> list[str]:
    """
    Returns a list of all references in the given
    databases.

    Parameters
    ----------
    dbs : Union[str, list[str]]
        One or more paths to SQLite databases produced by
        Nanocompore's eventalign_collapse or remora_resquiggle
        subcommands.
    has_data : bool, default=True
        If True (default) will return only
        references for which there are data.

    Returns
    -------
    list
        List of transcript reference id strings.

    Examples
    --------
    >>> from nanocompore.api import get_db_references
    >>> get_db_references('wt1.sqlite')
    ['ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', 'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|']
    """
    references = set()
    for db in dbs:
        if has_data:
            transcripts = PreprocessingDB(db).get_references_with_data()
            references.update(transcripts.keys())
        else:
            transcripts = PreprocessingDB(db).get_references()
            references.update([t.ref_id for t in transcripts])
    return list(references)


def _get_bam_reads(
        bams: Union[str, list[str]],
        reference_id: str,
        kit: Kit,
        selected_reads: Optional[list[str]]=None
    ) -> tuple[Float[np.ndarray, "reads positions variables"],
               list[str],
               list[str]]:
    """
    Get a numpy array with read data for specific reference_id from
    one or more BAMs produced by Uncalled4.

    Parameters
    ----------
    bams : str
        Paths to one or more BAM files.
    reference_id : str
        ID for a reference sequence (transcript).
    kit : Kit
        The chemistry used for the sequencing (e.g. RNA002 or RNA004).
    selected_reads : Optional[list[str]]
        Optional list of UUIDs of the reads for which to get data.
        By default it's set to None and returns all reads.

    Returns
    -------
    tuple[Float[np.ndarray, ["reads positions variables"]],
          list[str]]

        A tuple with (signal_data, reads, samples)

        - signal_data is a 3D array with shape (reads, positions, variables).
          In the variables dimension 0=intensity, 1=dwell time.
        - reads is a list of read ids (qname).
        - samples is a list of the sample ids (in the
          order in which the bams were provided).

    Raises
    ------
    KeyError
        If the reference_id is not found in the BAMs.

    Examples
    --------
    >>> from nanocompore.api import get_bam_reads, Kit
    >>> get_bam_reads(['path/to/wt.bam', 'path/to/kd.bam'], 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', Kit.RNA004)
    """
    if isinstance(bams, str):
        bams = [bams]
    opened_bams = {i: pysam.AlignmentFile(bam, 'rb')
                   for i, bam in enumerate(bams)}
    ref_len = None
    for bam in opened_bams.values():
        try:
            ref_len = bam.get_reference_length(reference_id)
            break
        except KeyError:
            pass
    else:
        raise KeyError(f"Reference {reference_id} not found in the bams.")

    uncalled4 = Uncalled4(reference_id,
                          ref_len,
                          opened_bams,
                          kit)
    data, reads, samples = uncalled4.get_data()
    if selected_reads:
        read_set = set(selected_reads)
        mask = np.array([read in read_set for read in reads])
        return data[mask], reads[mask].tolist(), samples[mask].tolist()
    return np.moveaxis(data, (0, 1, 2), (1, 0, 2)), reads.tolist(), samples.tolist()


def _get_db_reads(
        dbs: Union[str, list[str]],
        reference_id: str,
        selected_reads: Optional[list[str]]=None
    ) -> tuple[
            Float[np.ndarray, "reads positions variables"],
            list[str],
            list[str]]:
    """
    Get a numpy array with read data for specific reference_id from
    one or more preprocessing databases produced by eventalign_collapse or
    remora_resquiggle.

    Parameters
    ----------
    dbs : Union[str, list[str]]
        Paths to the SQLite database files.
    reference_id : str
        ID for a reference sequence (transcript).
    selected_reads : Optional[list[str]]
        Optional list of UUIDs of the reads for which to get data.
        By default it's set to None and returns all reads.

    Returns
    -------
    tuple[Float[np.ndarray, ["reads positions variables"]],
          list[str]]

        A tuple with (signal_data, reads, samples)

        - signal_data is a 3D array with shape (reads, positions, variables).
          In the variables dimension 0=intensity, 1=dwell time.
        - reads is a list of read ids (qname).
        - samples is a list of the sample ids (in the
          order in which the dbs were provided).

    Examples
    --------
    >>> from nanocompore.api import get_db_reads
    >>> get_db_reads(['path/to/wt.sqlite', 'path/to/kd.sqlite'], 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')
    """
    all_data = []
    all_reads = []
    all_samples = []
    for i, db in enumerate(dbs):
        data, reads = _get_sample_db_reads(db, reference_id, selected_reads)
        if len(reads) > 0:
            all_data.append(data)
            all_reads.extend(reads)
            all_samples.extend(np.full((len(reads),), i))
    return np.concatenate(all_data), all_reads, all_samples


def _get_sample_db_reads(
        db: str,
        reference_id: str,
        selected_reads: Optional[list[str]]=None
    ) -> tuple[
            Union[Float[np.ndarray, "reads positions variables"], None],
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
    tuple[Union[Float[np.ndarray, ["reads positions variables"]], None],
          list[str]]

        A tuple with (signal_data, reads)

        - signal_data is a 3D array with shape (reads, positions, variables).
          In the variables dimension 0=intensity, 1=dwell time. If no reads
          are found for the reference, signal_data will be None.
        - reads is a list of read ids (bam qname).

    Examples
    --------
    >>> from nanocompore.api import get_sample_db_reads
    >>> get_sample_db_reads('path/to/wt.sqlite', 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|')
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
        if len(reads) == 0:
            return None, reads
        signal_data = np.dstack([intensity, dwell])
        if selected_reads:
            read_set = set(selected_reads)
            mask = np.array([read in read_set for read in reads])
            filtered_reads = [read for read in reads if read in read_set]
            return signal_data[mask], filtered_reads
        return signal_data, reads


def _get_db_read(db: str, read_id: str) -> Float[np.ndarray, "reads positions variables"]:
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

    Examples
    --------
    >>> from nanocompore.api import get_db_read
    >>> get_db_read('path/to/wt.sqlite', 'f9733448-6e6b-47ba-9501-01eda2f5ea26')
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


def _get_bam_pos(bams: Union[str, list[str]], reference_id: str, pos: int, kit: Kit) -> pd.DataFrame:
    """
    Get data for a position from one or more BAMs created by Uncalled4.

    Returns the signal data for a specific position
    of the given reference transcript from all reads.

    Parameters
    ----------
    bams : Union[str, list[str]]
        Paths to the BAM files produced by Uncalled4.
    reference_id : str
        ID for a reference sequence (transcript).
    pos : int
        Position on the transcript for which to get data. A 0-based
        index is assumed.
    kit : Kit
        The chemistry used for the sequencing (e.g. RNA002 or RNA004).


    Returns
    -------
    pandas.DataFrame
        Where the DataFrame contains the following columns:

        - sample     id of the sample (equal to the index
                     of the bam in the bams parameter).
        - read:      id of the read
        - intensity: current intensity
        - dwell:     dwell time for the kmer

    Examples
    --------
    >>> from nanocompore.api import get_bam_pos
    >>> get_bam_pos(['wt1.bam', 'kd1.bam'], 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|', 532, Kit.RNA004)
       sample                                  read  intensity  dwell
    0       0  a4395b0d-dd3b-48e3-8afb-4085374b1147     3800.0    7.0
    1       0  f9733448-6e6b-47ba-9501-01eda2f5ea26     4865.0  126.0
    2       0  6f5e3b2e-f27b-47ef-b3c6-2ab4fdefd20a     3272.0   42.0
    3       1  3f46f499-8ce4-4817-8177-8ad61b784f27     4807.0   57.0
    4       1  73d62df4-f04a-4207-a4bc-7b9739b3c3b2     4336.0  132.0
    5       1  b7bc9a36-318e-4be2-a90f-74a5aa6439bf     -861.0    7.0
    """
    reads_data, reads, samples = _get_bam_reads(bams, reference_id, kit)
    pos_data = reads_data[:, pos, :]
    valid = ~np.isnan(pos_data[:, 0])
    return pd.DataFrame({'sample': np.array(samples)[valid],
                         'read': np.array(reads)[valid],
                         'intensity': pos_data[valid, INTENSITY_POS],
                         'dwell': pos_data[valid, DWELL_POS]})


def _get_db_pos(dbs: Union[str, list[str]], reference_id: str, pos: int) -> pd.DataFrame:
    """
    Get data for a position for a single sample from SQLite databases
    created by the eventalign_collapse or remora_resquiggle subcommands.

    Returns the signal data for a specific position
    of the given reference transcript from all reads.

    Parameters
    ----------
    dbs : Union[str, list[str]]
        Paths to the SQLite databases produced by eventalign_collapse
        or remora_resquiggle.
    reference_id : str
        ID for a reference sequence (transcript).
    pos : int
        Position on the transcript for which to get data. A 0-based
        index is assumed.

    Returns
    -------
    pandas.DataFrame
        Where the DataFrame contains the following columns:

        - sample     id of the sample (equal to the index
                     of the db in the dbs parameter).
        - read:      id of the read
        - intensity: current intensity
        - dwell:     dwell time for the kmer
    """
    reads_data, reads, samples = _get_db_reads(dbs, reference_id)
    pos_data = reads_data[:, pos, :]
    valid = ~np.isnan(pos_data[:, 0])
    return pd.DataFrame({'sample': np.array(samples)[valid],
                         'read': np.array(reads)[valid],
                         'intensity': pos_data[valid, INTENSITY_POS],
                         'dwell': pos_data[valid, DWELL_POS]})


def get_metadata(db: str) -> dict[str, str]:
    """
    Returns the metadata from the given SQLite database.

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


def _get_data_files(config: Config) -> dict[str, str]:
    field = 'db'
    if config.get_resquiggler() == UNCALLED4:
        field = 'bam'
    return {sample: sample_data[field]
            for cond_data in config.get_data().values()
            for sample, sample_data in cond_data.items()}

