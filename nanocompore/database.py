import os
import sqlite3
import json

from contextlib import closing
from pathlib import Path
from typing import Iterator, Union

import pandas as pd
import numpy as np

from loguru import logger

from nanocompore.common import NanocomporeError
from nanocompore.common import TranscriptRow


# ======= RESULT DATABASE QUERIES =======

CREATE_TRANSCRIPTS_RESULTS_TABLE = """
CREATE TABLE IF NOT EXISTS transcripts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name VARCHAR NOT NULL UNIQUE
);
"""
INSERT_TRANSCRIPTS_QUERY = """
INSERT INTO transcripts (id, name) VALUES(?, ?);
"""

CREATE_KMER_RESULTS_TABLE = """
CREATE TABLE IF NOT EXISTS kmer_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    transcript_id INTEGER NOT NULL,
    pos INTEGER NOT NULL,
    kmer INTEGER NOT NULL,
    UNIQUE (transcript_id, pos),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(id)
);
"""

CREATE_KMER_RESULTS_TRANSCRIPT_ID_INDEX = """
CREATE INDEX IF NOT EXISTS kmer_results_transcript_id_index
    ON kmer_results(transcript_id);
"""

CREATE_TRANSCRIPTS_NAME_INDEX = """
CREATE INDEX IF NOT EXISTS transcripts_name_index
    ON transcripts(name);
"""

# ======= PREPROCESSING DATABASES QUERIES =======

DROP_METADATA_TABLE_QUERY = """
DROP TABLE IF EXISTS metadata;
"""

CREATE_METADATA_TABLE_QUERY = """
CREATE TABLE metadata (
    key TEXT NOT NULL,
    value TEXT NOT NULL,
    PRIMARY KEY(key)
);
"""

DROP_SIGNAL_DATA_TABLE_QUERY = """
DROP TABLE IF EXISTS signal_data;
"""

CREATE_SIGNAL_DATA_TABLE_QUERY = """
CREATE TABLE signal_data (
    transcript_id INTEGER NOT NULL,
    read_id INTEGER NOT NULL,
    intensity BLOB,
    dwell BLOB
);
"""

INSERT_SIGNAL_DATA_QUERY = """
INSERT INTO signal_data VALUES (?, ?, ?, ?);
"""

DROP_READS_TABLE_QUERY = """
DROP TABLE IF EXISTS reads;
"""

CREATE_READS_TABLE_QUERY = """
CREATE TABLE reads (
    read TEXT NOT NULL,
    id INTEGER NOT NULL,
    invalid_kmers REAL
);
"""

DROP_READS_ID_INDEX_QUERY = """
DROP INDEX IF EXISTS reads_id_index;
"""

CREATE_READS_ID_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS reads_id_index ON reads (id);
"""

DROP_SIGNAL_DATA_INDEX_QUERY = """
DROP INDEX IF EXISTS signal_data_index;
"""

CREATE_SIGNAL_DATA_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS signal_data_index ON signal_data (transcript_id, read_id);
"""

DROP_TRANSCRIPTS_TABLE_QUERY = """
DROP TABLE IF EXISTS transcripts;
"""

CREATE_TRANSCRIPTS_TABLE_QUERY = """
CREATE TABLE transcripts (
    name TEXT NOT NULL,
    id INTEGER NOT NULL,
    PRIMARY KEY(name)
);
"""

DROP_TRANSCRIPTS_ID_INDEX_QUERY = """
DROP INDEX IF EXISTS transcripts_id_index;
"""

CREATE_TRANSCRIPTS_ID_INDEX_QUERY = """
CREATE INDEX IF NOT EXISTS transcripts_id_index ON transcripts (id);
"""

INSERT_READS_QUERY = """
INSERT INTO reads (read, id, invalid_kmers) VALUES(?, ?, ?);
"""

GET_SIGNAL_DATA_FOR_TRANSCRIPT_QUERY = """
SELECT intensity, dwell
FROM signal_data sd
JOIN reads r ON sd.read_id = r.id
JOIN transcripts t ON sd.transcript_id = t.id
WHERE t.name = ?
  AND r.invalid_kmers <= ?
"""

BASE_KMER_RESULT_COLUMNS = ['id', 'transcript_id', 'pos', 'kmer']


class ResultsDB():
    def __init__ (self, config, init_db=False):
        self._config = config

        if config.get_outpath():
            self._outpath = config.get_outpath()
        else:
            self._outpath = os.getcwd()

        self._prefix = config.get_outprefix()
        self.db_path = os.path.join(self._outpath, f"{self._prefix}sampComp_sql.db")
        if init_db:
            self._setup_database(config.get_result_exists_strategy())
            self._create_tables()


    def get_results(self,
                    columns: Union[list[str], None]=None,
                    chunksize: Union[int, None]=None) -> Union[pd.DataFrame, Iterator[pd.DataFrame]]:
        if not columns:
            cols = '*'
        else:
            cols = ','.join(columns)
        with closing(sqlite3.connect(self.db_path)) as conn:
            query = f"""
            SELECT {cols}
            FROM kmer_results res
            JOIN transcripts t ON t.id = res.transcript_id
            """
            for chunk in pd.read_sql(query, conn, chunksize=chunksize):
                yield chunk


    def get_all_results(self):
        with closing(sqlite3.connect(self.db_path)) as conn:
            query = """
            SELECT *
            FROM kmer_results res
            JOIN transcripts t ON t.id = res.transcript_id
            """
            return pd.read_sql(query, conn)


    def get_transcripts(self):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            query = "SELECT name FROM transcripts"
            return [row[0] for row in cursor.execute(query).fetchall()]


    def get_result_column_names(self) -> list[str]:
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            res = cursor.execute('SELECT * FROM kmer_results LIMIT 1')
            return [description[0] for description in res.description]


    def get_columns(self, cols: list[str]) -> pd.DataFrame:
        col_list = ', '.join(cols)
        with closing(sqlite3.connect(self.db_path)) as conn:
            query = f"""
            SELECT {col_list}
            FROM kmer_results res
            JOIN transcripts t ON t.id = res.transcript_id
            """
            return pd.read_sql(query, conn)


    def get_columns_for_ref(self, cols: list[str], ref: str) -> pd.DataFrame:
        """
        Returns the data for the given reference, selecting
        only the listed columns.

        Parameters
        ----------
        cols : list[str]
            List of column names to select.
        ref : str
            Reference id.

        Returns
        -------
        pd.DataFrame
            The dataframe would contain a "pos" column
            along with the list of provided columns.
        """
        cols = ['pos'] + cols
        col_list = ', '.join(cols)
        query = f'''
        SELECT {col_list}
        FROM kmer_results res
        JOIN transcripts t ON t.id = res.transcript_id
        WHERE t.name = ?
        '''
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            results = cursor.execute(query, (ref,)).fetchall()

            df = pd.DataFrame([dict(zip(cols, row))
                               for row in results])
            return df


    def get_next_transcript_id(self):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            query = "SELECT MAX(id) FROM transcripts"
            return cursor.execute(query).fetchone()[0] + 1


    def save_test_results(self, transcript, test_results):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(INSERT_TRANSCRIPTS_QUERY, (transcript.id, transcript.name))
            test_columns = dict(zip(test_results.columns, test_results.dtypes))
            self._create_missing_columns(test_columns, cursor)
            test_results.to_sql('kmer_results', conn, if_exists='append', index=False)


    def save_transcript(self, transcript):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(INSERT_TRANSCRIPTS_QUERY, (transcript.id, transcript.name))
            conn.commit()


    def index_database(self):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_TRANSCRIPTS_NAME_INDEX)
            cursor.execute(CREATE_KMER_RESULTS_TRANSCRIPT_ID_INDEX)


    def _create_missing_columns(self, test_columns, cursor):
        info = cursor.execute("SELECT * FROM kmer_results LIMIT 1")
        existing_columns = {item[0] for item in info.description}
        logger.trace("Adding columns to table kmer_results")
        try:
            for column, column_type in test_columns.items():
                if column in existing_columns:
                    continue

                if column_type in (float, np.float64, np.float32):
                    column_type = 'FLOAT'
                elif column_type is str:
                    column_type = 'VARCHR'
                elif column_type is int:
                    column_type = 'INTEGER'
                else:
                    column_type = 'VARCHAR'

                update_query = f"ALTER TABLE kmer_results ADD COLUMN {column} {column_type}"
                cursor.execute(update_query)
                existing_columns.add(column)
                logger.trace(f"Added {column} to kmer_results")
        except Exception:
            raise NanocomporeError(f"Database error: Failed to insert at least one of {test_columns} new labels into kmer_results of {self.db_path}")


    def _setup_database(self, result_exists_strategy):
        if os.path.isfile(self.db_path):
            if result_exists_strategy == 'overwrite':
                os.remove(self.db_path)
                logger.debug(f"Removed existing database file '{self.db_path}'")
            elif result_exists_strategy == 'continue':
                logger.info(f"Database file '{self.db_path}' already exists and result_exists_strategy is set to 'continue'. Will try to reuse it.")
            else:
                raise NanocomporeError(f"Database file '{self.db_path}' exists and 'results_exists_strategy' is 'stop'")
        with closing(sqlite3.connect(self.db_path)) as conn:
            conn.execute('PRAGMA foreign_keys = ON')
            conn.execute('PRAGMA journal_mode = wal')
            conn.execute('PRAGMA synchronous = NORMAL')


    def _create_tables(self):
        with closing(sqlite3.connect(self.db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_TRANSCRIPTS_RESULTS_TABLE)
            cursor.execute(CREATE_KMER_RESULTS_TABLE)


class PreprocessingDB:
    def __init__(self, db):
        self._db = db


    def connect(self):
        conn = sqlite3.connect(self._db, isolation_level=None)
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = OFF')
        # Use largest possible page size of 65kb
        conn.execute('PRAGMA page_size = 65536')
        # Use cache of 2048 x page_size. That's 128MB of cache.
        conn.execute('PRAGMA cache_size = 2000')
        return conn


    def setup(self, metadata):
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(DROP_READS_ID_INDEX_QUERY)
            cursor.execute(DROP_SIGNAL_DATA_INDEX_QUERY)
            cursor.execute(DROP_SIGNAL_DATA_TABLE_QUERY)
            cursor.execute(CREATE_SIGNAL_DATA_TABLE_QUERY)
            cursor.execute(DROP_READS_TABLE_QUERY)
            cursor.execute(CREATE_READS_TABLE_QUERY)
            cursor.execute(DROP_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(DROP_METADATA_TABLE_QUERY)
            cursor.execute(CREATE_METADATA_TABLE_QUERY)
            cursor.execute(DROP_TRANSCRIPTS_ID_INDEX_QUERY)
        self.write_metadata(metadata)


    def create_indices(self):
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_READS_ID_INDEX_QUERY)
            cursor.execute(CREATE_SIGNAL_DATA_INDEX_QUERY)
            cursor.execute(CREATE_TRANSCRIPTS_ID_INDEX_QUERY)


    def write_metadata(self, metadata):
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            for key, value in metadata.items():
                cursor.execute("INSERT INTO metadata (key, value) VALUES (?, ?)",
                               (key, str(value)))


    def merge_in_databases(self, databases, merge_reads=True):
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            for other_db in databases:
                cursor.execute(f"ATTACH '{other_db}' as tmp")
                cursor.execute("BEGIN")
                cursor.execute("INSERT INTO transcripts SELECT * FROM tmp.transcripts")
                if merge_reads:
                    cursor.execute("INSERT INTO reads SELECT * FROM tmp.reads")
                cursor.execute("INSERT INTO signal_data SELECT * FROM tmp.signal_data")
                conn.commit()
                cursor.execute("DETACH DATABASE tmp")
                # Delete the other db
                Path(other_db).unlink()


    def get_references_with_data(self) -> dict[str, int]:
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            query = """
            SELECT t.name, COUNT(sd.read_id)
            FROM signal_data sd
            INNER JOIN transcripts t ON sd.transcript_id = t.id
            GROUP BY t.id
            """
            return {row[0]: row[1]
                    for row in cursor.execute(query).fetchall()}


    def get_references(self):
        with closing(self.connect()) as conn,\
             closing(conn.cursor()) as cursor:
            query = """
            SELECT DISTINCT name, id
            FROM transcripts
            """
            return {TranscriptRow(row[0], row[1])
                    for row in cursor.execute(query).fetchall()}


    @staticmethod
    def write_signal_data_rows(connection, rows):
        with closing(connection.cursor()) as cursor:
            cursor.execute("begin")
            cursor.executemany(INSERT_SIGNAL_DATA_QUERY, rows)
            cursor.execute("commit")


    @staticmethod
    def get_signal_data(connection,
                        transcript_name: str,
                        max_invalid_ratio: float):
        """
        Get signal data for a transcript.
        Will return reads with lower invalid ratio first.

        Parameters
        ----------
        connection : sqlite3.connection
            Connection to the database.
        transcript_name : str
            The name of the transcript (transcript reference).
        max_invalid_ratio : float
            Will return only reads with invalid ratio smaller
            than the given number.

        Returns
        -------
        list[tuple[bytearray, bytearray]]
            List of tuples (intensity, dwell), where
            intensity and dwell are binary arrays
            than need to be decoded to numpy arrays.
        """
        with closing(connection.cursor()) as cursor:
            return cursor.execute(GET_SIGNAL_DATA_FOR_TRANSCRIPT_QUERY,
                                  (transcript_name, max_invalid_ratio)).fetchall()

