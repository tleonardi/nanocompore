import sqlite3 as sql
import os, collections, sys

from loguru import logger
import pandas as pd
import numpy as np

from nanocompore.common import *


class SampComp_DB():

    def __init__ (self, outpath='', prefix='', result_exists_strategy="stop"):
        if outpath:
            self._outpath = outpath
        else:
            self._outpath = os.getcwd()

        self._prefix = prefix
        self._db_path = os.path.join(self._outpath, f"{self._prefix}sampComp_sql.db")

        self._enterDB(result_exists_strategy)

        self._create_tables()


    def _checkIfLabelsExist(self, test_results, table=''):
        if not table:
            logger.debug(f"{table} doesn't exist in the database yet. Creating")
            table = self._default_storage_table

        labels = self._create_test_labels(test_results)

        info = self._cursor.execute(f"SELECT * FROM {table}")
        columns = [item[0] for item in info.description]
        logger.trace(f"Adding data to table {table}")
        try:
            for label, label_type in labels:
                if label not in columns:
                    if label_type == float or label_type == np.float64:
                        label_type = 'FLOAT'
                    elif label_type == str:
                        label_type = 'VARCHR'
                    elif label_type == int:
                        label_type = 'INTEGER'
                    else:
                        label_type = 'VARCHAR'

                    update_query = f"ALTER TABLE {table} ADD COLUMN {label} {label_type}"
                    self._connection.execute(update_query)
                    columns.append(label)
                    self._connection.commit()
                    logger.trace(f"Added {label} to {table}")
        except:
            raise NanocomporeError(f"Database error: Failed to insert at least one of {labels} new labels into {table} of {self._db_path}")


    def _create_test_labels(self, test_results):
        labels = set()
        for pos in test_results:
            for test in test_results[pos]:
                test_type = type(test_results[pos][test])
                labels.add((test, test_type))
        logger.trace(f'tests and test labels created for the table')
        return labels


    def _makes_values(self, pos=0, tx_id=0, results_dict={}, kmer_table_headers=[]):
        if results_dict:
            logger.trace(f"Getting the kmer_seq id for each position of the ref_id")
            kmer_seq_id = self._insert_query_and_get_id(query=results_dict['kmer_seq'], table='kmer_seqs', column='sequence')
            logger.trace(f"The kmer {results_dict['kmer_seq']} has the kmer_id {kmer_seq_id}")
            values = [tx_id, pos, kmer_seq_id]
            for test_type in kmer_table_headers:
                if test_type in results_dict and results_dict[test_type]:
                    values.append(results_dict[test_type])
                else:
                    values.append(np.nan)
            return values
        else:
            raise NanocomporeError(f"No results for position {pos}")


    def store_test_results(self, tx_name, test_results, table=''):
        if self._is_connected():
            tx_id = self._insert_query_and_get_id(query=tx_name, table='transcripts', column='name')
            self._checkIfLabelsExist(test_results, table=table)

            kmer_table_headers = self.get_kmer_table_headers()
            if len(kmer_table_headers) < 1:
                raise NanocomporeError('DataBaseError: Error collecting kmer_stats column headers')
            else:
                table_cols = "transcript_id, pos, kmer_seq_id, {}".format(', '.join(kmer_table_headers))
            for pos in sorted(test_results):
                results_dict = test_results[pos]
                try:
                    values=self._makes_values(pos, tx_id, results_dict, kmer_table_headers)
                except:
                    raise NanocomporeError("DataBaseError: Error gathering values for entry into the database")

                try:
                    table_values = ', '.join(['?'] * len(values))
                    insert_querry = f"INSERT INTO kmer_stats ({table_cols}) VALUES ({table_values})"
                    self._cursor.execute(insert_querry, values)
                except:
                    self._exitDB()
                    raise NanocomporeError(f"Error storing statistics for transcript '{tx_name}', pos {pos}")
            self._connection.commit()
            logger.debug(f"Results for {tx_name} added to the database")
        else:
            raise NanocomporeError("Not connected to the Database")


    def store_read_level_data(self, tx_name, read_data):
        """
        Stores the measurements (intensity and dwell time) for
        each position of every read to the database.
        """
        if self._is_connected():
            # Make sure that the transcript exists in the database and get its id
            tx_id = self._insert_query_and_get_id(query=tx_name, table='transcripts', column='name')
            # Make sure that the kmers exist in the database and get their ids
            kmer_to_id = self._get_or_create_kmer_ids('kmer_seqs', 'sequence', read_data['kmer'])

            # Write the read ids to the reads table
            reads = read_data[['sample', 'read']].copy()
            reads['transcript_id'] = tx_id
            reads.rename(columns={'read': 'id'}, inplace=True)
            self._insert_missing('reads', reads.columns, list(reads.itertuples(index=False)))

            # Write the intensity/dwell time data for all read/position pairs
            # to the read_level_data table
            read_data['kmer_seq_id'] = [kmer_to_id[x] for x in read_data['kmer']]
            read_data.drop(columns=['kmer'], inplace=True)
            read_data.rename(columns={'read': 'read_id'}, inplace=True)
            read_data.to_sql('read_level_data',
                             self._connection,
                             if_exists='append',
                             index=False,
                             method='multi',
                             chunksize=32000)

            logger.debug(f"Read level data for {tx_name} added to the database")
        else:
            raise NanocomporeError("Not connected to the Database")


    def get_kmer_table_headers(self):
        table_headers = self._cursor.execute('''SELECT * FROM kmer_stats''')
        headers = []
        for x in table_headers.description:
            header = x[0]
            if header not in ['id', 'pos', 'transcript_id', 'kmer_seq_id']:
                headers.append(header)
        return headers


    def getAllData(self):
        data_query = self._make_query()
        data = pd.read_sql(data_query, self._connection)

        return data


    def closeDB(self):
        self._exitDB()


    def get_transcripts(self):
        return [row[0]
                for row in self._cursor.execute("SELECT name FROM transcripts").fetchall()]


    ################## Private methods ##################
    def _enterDB(self, result_exists_strategy):
        if result_exists_strategy == 'overwrite':
            if os.path.isfile(self._db_path):
                os.remove(self._db_path)
                #sys.stderr.write(f"Removed existing database file '{self._db_path}'\n")
                logger.debug(f"Removed existing database file '{self._db_path}'")
            self._connect_to_db()
        elif result_exists_strategy == 'continue':
            logger.info(f"Database file '{self._db_path}' already exists and result_exists_strategy is set to 'continue'. Will try to reuse it.")
            self._connect_to_db()
        else:
            if os.path.isfile(self._db_path):
                raise NanocomporeError(f"database file '{self._db_path}' exists and overwrite is False")
            else:
                self._connect_to_db()

    def _connect_to_db(self):
        self._connection = sql.connect(self._db_path)
        self._connection.execute('PRAGMA foreign_keys = ON')
        self._connection.execute('PRAGMA journal_mode = wal')
        self._connection.execute('PRAGMA synchronous = NORMAL')
        self._cursor = self._connection.cursor()
        logger.debug(f"Connected to {self._db_path}")

    def _exitDB(self):
        if self._is_connected():
            self._cursor.close()
            self._connection.close()
            logger.debug(f"Closed the connection to {self._db_path}")
        else:
            logger.debug(f"{self._db_path} is not connected")

    def _is_connected(self):
     try:
        self._connection.cursor()
        return True
     except Exception as ex:
        return False

    def _create_tables(self):
        # TODO
        # "parameters" table:
        '''
        table_def_parameters = ["univariate_test VARCHAR CHECK (univariate_test in ('ST', 'MW', 'KS'))",
                                "gmm_covariance_type VARCHAR",
                                "gmm_test VARCHAR CHECK (gmm_test in ('anova', 'logit'))"]
        '''

        # "transcripts" table:
        table_def_transcripts = ["id INTEGER PRIMARY KEY AUTOINCREMENT",
                                "name VARCHAR NOT NULL UNIQUE"]

        #TODO
        # "samples" table:
        '''
        table_def_samples = ["id INTEGER PRIMARY KEY AUTOINCREMENT",
                             "name VARCHAR NOT NULL UNIQUE",
                             "db_path VARCHAR NOT NULL"]
        '''

        #TODO
        # "kmer_seqs" table:
        table_def_kmer_seqs = ["id INTEGER PRIMARY KEY AUTOINCREMENT",
                               "sequence VARCHAR NOT NULL UNIQUE"]

        #TODO
        # "whitelist" table:
        '''
        table_def_whitelist = ["transcript_id INTEGER NOT NULL",
                            "sample_id INTEGER NOT NULL UNIQUE",
                            "read_id INTEGER NOT NULL", # foreign key for "reads" table in EventAlign DB
                            "UNIQUE (sample_id, read_id)"
                            "FOREIGN KEY (transcript_id) REFERENCES transcripts(id)",
                            "FOREIGN KEY (sample_id) REFERENCES samples(id)"]
        '''

        # "kmer_stats" table:
        table_def_kmer_stats = ["id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL",
                                "transcript_id INTEGER NOT NULL",
                                "kmer_seq_id INTEGER NOT NULL",
                                "pos INTEGER NOT NULL",
                                "UNIQUE (transcript_id, pos)",
                                "FOREIGN KEY(transcript_id) REFERENCES transcripts(id)",
                                "FOREIGN KEY(kmer_seq_id) REFERENCES kmer_seqs(id)"]

        # "reads" table:
        table_def_reads = ["sample VARCHAR NOT NULL",
                           "id VARCHAR NOT NULL",
                           "transcript_id INTEGER NOT NULL",
                           "PRIMARY KEY (sample, id)"]
        # TODO: maybe add metadata for the read (e.g. quality, cigar string, etc.)

        table_read_level_data = [
            'condition VARCHAR NOT NULL',
            'sample VARCHAR NOT NULL',
            'read_id INTEGER NOT NULL',
            'pos INTEGER NOT NULL',
            'kmer_seq_id INTEGER NOT NULL',
            'intensity FLOAT NOT NULL',
            'dwell FLOAT NOT NULL',
            'PRIMARY KEY (condition, sample, read_id, pos)',
            'FOREIGN KEY(sample, read_id) REFERENCES reads(sample, id)',
            'FOREIGN KEY(kmer_seq_id) REFERENCES kmer_seqs(id)']

        table_defs = collections.OrderedDict([('transcripts', table_def_transcripts),
                                              ('kmer_seqs', table_def_kmer_seqs),
                                              ('kmer_stats', table_def_kmer_stats),
                                              ('reads', table_def_reads),
                                              ('read_level_data', table_read_level_data)])

        for table in table_defs:
            table_headers = ', '.join(table_defs[table])
            construction_querry = f"CREATE TABLE IF NOT EXISTS {table} ({table_headers})"
            error = None
            for retry in range(3):
                try:
                    self._cursor.execute(construction_querry)
                    self._connection.commit()
                    logger.debug(f"created {table}")
                    break
                except Exception as e:
                    logger.error(f"Error creating {table} table. Retrying.")
                    error = e
            else:
                logger.error(f"Error creating {table} table after 3 retries.")
                raise error

        self._default_storage_table = 'kmer_stats'


    def _get_or_create_kmer_ids(self, table, column, values):
        """
        Gets the ids of the values in the table or creates them if they don't exist.
        Works only only for tables with a single column."""
        insert_query =  f"INSERT INTO {table} ({column}) VALUES (?) ON CONFLICT DO NOTHING"
        self._cursor.executemany(insert_query, [(x, ) for x in values])
        seq_to_id = {}
        # The SQLite driver supports maximum 999 parameters
        # so we need to split the query into batches
        for i in range(0, len(values), 999):
            batch = values[i:i+999]
            param_definition = ', '.join(['?'] * len(batch))
            select_query = f"SELECT id, {column} FROM {table} WHERE {column} IN ({param_definition})"
            new_seqs = self._cursor.execute(select_query, tuple(batch)).fetchall()
            seq_to_id.update({x[1]: x[0] for x in new_seqs})
        return seq_to_id


    def _insert_missing(self, table, columns, values):
        """
        Inserts the values into the table if they don't exist.
        """
        insert_query =  f"""
        INSERT INTO {table} ({','.join(columns)})
        VALUES ({','.join(['?'] * len(columns))})
        ON CONFLICT DO NOTHING
        """
        self._cursor.executemany(insert_query, values)


    def _insert_query_and_get_id(self, query='', table='', column=''):
        if all(map(lambda x: bool(x), [query, table, column])):
            error = None
            for retry in range(3):
                try:
                    self._cursor.execute(f"SELECT id FROM {table} WHERE {column} = ?", [query])
                    row = self._cursor.fetchone()
                    if row is not None:
                        return row[0]
                    else:
                        self._cursor.execute(f"INSERT INTO {table} ({column}) VALUES (?)", [query])
                        self._connection.commit()
                        return self._cursor.lastrowid
                except Exception as e:
                    logger.error(f"Failed to insert/look up {table}, {column}, {query}. Retrying.")
                    error = e
            else:
                raise NanocomporeError(f"Failed to insert/look up {table}, {column}, {query} after 3 retries.")
        else:
            raise NanocomporeError(f"At least one input was empty query={query}, table={table}, column={column}")


    def _make_query(self):
        table_cols = ['kmer_stats.pos', 'tx.name', 'kmer_seqs.sequence']

        for header in self.get_kmer_table_headers():
            table_cols.append(f"kmer_stats.{header}")
        table_cols = ', '.join(table_cols)

        data_query = f"SELECT {table_cols} FROM kmer_stats kmer_stats LEFT JOIN transcripts tx ON tx.id = kmer_stats.transcript_id LEFT JOIN kmer_seqs kmer_seqs ON kmer_stats.kmer_seq_id = kmer_seqs.id"
        return data_query

    def _index_database(self):
        index_queries = []
        index_queries.append("CREATE INDEX IF NOT EXISTS kmer_stats_tx_id_index ON kmer_stats (transcript_id)")
        index_queries.append("CREATE INDEX IF NOT EXISTS kmer_stats_kmer_id_index ON kmer_stats (kmer_seq_id)")
        index_queries.append("CREATE UNIQUE INDEX IF NOT EXISTS transcripts_index ON transcripts (id, name)")
        index_queries.append("CREATE UNIQUE INDEX IF NOT EXISTS kmer_seqs_index ON kmer_seqs (id, sequence)")
        index_queries.append("CREATE INDEX IF NOT EXISTS read_level_data_read_id_index ON read_level_data (read_id)")
        index_queries.append("CREATE INDEX IF NOT EXISTS reads_transcript_id_index ON reads (transcript_id)")

        for index_query in index_queries:
            try:
                self._cursor.execute(index_query)
            except:
                index = index_query.split('ON')[0].split('INDEX')[1].strip()
                raise NanocomporeError(f"Error creating index {index}")
