import os
import sqlite3

from contextlib import closing

import pandas as pd
import numpy as np

from loguru import logger

from nanocompore.common import encode_kmer
from nanocompore.common import NanocomporeError


CREATE_TRANSCRIPTS_TABLE = """
CREATE TABLE IF NOT EXISTS transcripts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name VARCHAR NOT NULL UNIQUE
);
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

BASE_KMER_RESULT_COLUMNS = ['id', 'transcript_id', 'pos', 'kmer']


class SampCompDB():
    def __init__ (self, outpath, prefix, result_exists_strategy):
        if outpath:
            self._outpath = outpath
        else:
            self._outpath = os.getcwd()

        self._prefix = prefix
        self._db_path = os.path.join(self._outpath, f"{self._prefix}sampComp_sql.db")
        self._setup_database(result_exists_strategy)
        self._create_tables()


    def get_all_results(self):
        with closing(sqlite3.connect(self._db_path)) as conn:
            query = """
            SELECT *
            FROM kmer_results res
            JOIN transcripts t ON t.id = res.transcript_id
            """
            return pd.read_sql(query, conn)


    def get_transcripts(self):
        with closing(sqlite3.connect(self._db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            query = "SELECT name FROM transcripts"
            return [row[0] for row in cursor.execute(query).fetchall()]


    def save_test_results(self, transcript, test_results):
        with closing(sqlite3.connect(self._db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute("INSERT INTO transcripts (id, name) VALUES (?, ?)", (transcript.id, transcript.name))
            # test_columns = self._get_test_columns(test_results)
            test_columns = dict(zip(test_results.columns, test_results.dtypes))
            # self._create_missing_columns(test_columns, cursor)
            self._create_missing_columns(test_columns, cursor)
            test_results.to_sql('kmer_results', conn, if_exists='append', index=False)


            # query, data = self._prepare_test_results_query_and_data(tx_id,
            #                                                         test_results,
            #                                                         test_columns)
            # cursor.execute('commit')
            # cursor.execute('begin')
            # cursor.executemany(query, data)
            # cursor.execute('commit')


    def index_database(self):
        with closing(sqlite3.connect(self._db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_TRANSCRIPTS_NAME_INDEX)
            cursor.execute(CREATE_KMER_RESULTS_TRANSCRIPT_ID_INDEX)


    def _create_missing_columns(self, test_columns, cursor):
        info = cursor.execute(f"SELECT * FROM kmer_results LIMIT 1")
        existing_columns = {item[0] for item in info.description}
        logger.trace(f"Adding columns to table kmer_results")
        try:
            for column, column_type in test_columns.items():
                if column in existing_columns:
                    continue

                if column_type == float or column_type == np.float64:
                    column_type = 'FLOAT'
                elif column_type == str:
                    column_type = 'VARCHR'
                elif column_type == int:
                    column_type = 'INTEGER'
                else:
                    column_type = 'VARCHAR'

                update_query = f"ALTER TABLE kmer_results ADD COLUMN {column} {column_type}"
                cursor.execute(update_query)
                existing_columns.add(column)
                logger.trace(f"Added {column} to kmer_results")
        except:
            raise NanocomporeError(f"Database error: Failed to insert at least one of {test_columns} new labels into kmer_results of {self._db_path}")


    def _get_test_columns(self, test_results):
        return {key: type(value)
                for pos_results in test_results.values()
                for key, value in pos_results.items()}


    # def _prepare_test_results_query_and_data(self, tx_id, test_results, test_columns):
    #     test_columns = list(test_columns.keys())
    #     columns = ['transcript_id', 'pos', 'kmer'] + test_columns
    #     params = ['?' for _ in columns]
    #     data = [(tx_id,
    #              pos,
    #              encode_kmer(pos_results['kmer']),
    #              *[pos_results.get(col, None) for col in test_columns])
    #             for pos, pos_results in test_results.items()]
    #     query = f"""
    #     INSERT INTO kmer_results ({','.join(columns)})
    #     VALUES ({','.join(params)})"""
    #     return query, data


    def _setup_database(self, result_exists_strategy):
        if os.path.isfile(self._db_path):
            if result_exists_strategy == 'overwrite':
                os.remove(self._db_path)
                logger.debug(f"Removed existing database file '{self._db_path}'")
            elif result_exists_strategy == 'continue':
                logger.info(f"Database file '{self._db_path}' already exists and result_exists_strategy is set to 'continue'. Will try to reuse it.")
            else:
                raise NanocomporeError(f"Database file '{self._db_path}' exists and 'results_exists_strategy' is 'stop'")
        with closing(sqlite3.connect(self._db_path)) as conn:
            conn.execute('PRAGMA foreign_keys = ON')
            conn.execute('PRAGMA journal_mode = wal')
            conn.execute('PRAGMA synchronous = NORMAL')


    def _create_tables(self):
        with closing(sqlite3.connect(self._db_path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_TRANSCRIPTS_TABLE)
            cursor.execute(CREATE_KMER_RESULTS_TABLE)

