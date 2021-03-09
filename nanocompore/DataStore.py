# -*- coding: utf-8 -*-

from collections import *
import datetime
import os
import sqlite3 as lite

# Third party
from loguru import logger
import nanocompore as pkg

class DataStore(object):
    """ Init analysis and check args"""

    create_reads_query = ("CREATE TABLE IF NOT EXISTS reads ("
                          "id INTEGER NOT NULL PRIMARY KEY,"
                          "name VARCHAR NOT NULL UNIQUE,"
                          "sampleid INTEGER NOT NULL,"
                          "transcriptid VARCHAR NOT NULL,"
                          "refstart INT NOT NULL,"
                          "refend INT NOT NULL,"
                          "numevents INT NOT NULL,"
                          "numsignals INT NOT NULL,"
                          "dwelltime REAL NOT NULL,"
                          "FOREIGN KEY(sampleid) REFERENCES samples(id)"
                          "FOREIGN KEY(transcriptid) REFERENCES transcripts(id),"
                          "UNIQUE(id, name)"
                          ")"
                          )

    create_kmers_query = ("CREATE TABLE IF NOT EXISTS kmers ("
                          "id INTEGER NOT NULL PRIMARY KEY,"
                          "readid INTEGER NOT NULL,"
                          "position INTEGER NOT NULL,"
                          "sequence INTEGER NOT NULL,"
                          "num_events INTEGER NOT NULL,"
                          "num_signals INTEGER NOT NULL,"
                          "dwell_time REAL NOT NULL,"
                          "NNNNN_dwell_time REAL NOT NULL,"
                          "mismatch_dwell_time REAL NOT NULL,"
                          "median REAL NOT NULL,"
                          "mad REAL NOT NULL,"
                          "FOREIGN KEY(readid) REFERENCES reads(id)"
                          ")"
                          )
    # TODO: 'sequence' is stored redundantly - move it to a separate table


    create_samples_query = ("CREATE TABLE IF NOT EXISTS samples ("
                            "id INTEGER NOT NULL PRIMARY KEY,"
                            "name VARCHAR NOT NULL UNIQUE"
                            ")"
                            )


    create_transcripts_query = ("CREATE TABLE IF NOT EXISTS transcripts ("
                                "id INTEGER NOT NULL PRIMARY KEY,"
                                "name VARCHAR NOT NULL UNIQUE"
                                ")"
                                )


    def __init__(self, db_path:str):
        self.__db_path=db_path
        db_is_new = not os.path.exists(self.__db_path)
        logger.debug(f"DB file doesn't exist: {db_is_new}")
        if db_is_new: self.__init_db()

    def __enter__(self):
        self.__open_db_connection()
        return self

    def __exit__(self,exc_type, exc_value, traceback):
        self.__connection.commit()
        self.__close_db_connection()

    def __open_db_connection(self):
        try:
            logger.debug("Connecting to DB")
            self.__connection = lite.connect(self.__db_path);
            self.__cursor = self.__connection.cursor()
        except:
            logger.error("Error connecting to database")
            raise

    def __close_db_connection(self):
        if self.__connection:
            logger.debug("Closing connection to DB")
            self.__connection.commit()
            self.__connection.close()
            self.__connection = None
            self.__cursor = None

    def __init_db(self):
        logger.debug("Setting up DB tables")
        self.__open_db_connection()
        self.__cursor.execute(self.create_reads_query)
        self.__cursor.execute(self.create_kmers_query)
        self.__cursor.execute(self.create_samples_query)
        self.__cursor.execute(self.create_transcripts_query)
        self.__connection.commit()
        self.__close_db_connection()

    def store_read(self, read):
        """
        Insert data corresponding to a read into the DB.
        Args:
            read (Read): an instance of class Read that contains read-level data.
        Returns:
            Bool: returns True is read added successfully.
        """
        tx_id = self.get_transcript_id_by_name(read.ref_id, create_if_not_exists=True)
        sample_id = self.get_sample_id_by_name(read.sample_name, create_if_not_exists=True)
        self.__cursor.execute("INSERT INTO reads VALUES(NULL, ?, ?, ?, ?, ?, ?, ?, ?)",
                              (read.read_id, sample_id, tx_id, read.ref_start, read.ref_end,
                               read.n_events, read.n_signals, read.dwell_time))
        read_id = self.__cursor.lastrowid
        for kmer in read.kmer_l:
            self.__store_kmer(kmer=kmer, read_id=read_id)
        # TODO check for success and return true/false


    def __store_kmer(self, kmer, read_id):
        """
        Insert data corresponding to a kmer into the DB.
        Args:
            kmer (Kmer): instance of class Kmer that contains kmer-level data
            read_id (int): DB key of the read the kmer belongs to
        """
        res = kmer.get_results() # needed for 'median' and 'mad' values
        self.__cursor.execute("INSERT INTO kmers VALUES(NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                              (read_id, res["ref_pos"], res["ref_kmer"], res["num_events"],
                               res["num_signals"], res["dwell_time"], res["NNNNN_dwell_time"],
                               res["mismatch_dwell_time"], res["median"], res["mad"]))
        # TODO: should 'res["status"]' get stored as well?
        # TODO: check for success?


    def get_transcript_id_by_name(self, tx_name, create_if_not_exists=False):
        # TODO: This function should cache results
        if create_if_not_exists:
            query = ("INSERT INTO transcripts(id, name) "
                     f"SELECT NULL,'{tx_name}' "
                     " WHERE NOT EXISTS ( "
                     "  SELECT 1"
                     "  FROM transcripts"
                     f"  WHERE name = '{tx_name}' "
                     ");"
                     )
            try:
                self.__cursor.execute(query)
            except Exception:
                logger.error("There was an error while inserting a new transcript in the DB")
                raise Exception

        query = f"SELECT id from transcripts WHERE name = '{tx_name}'"
        try:
            self.__cursor.execute(query)
            record = self.__cursor.fetchone()
            self.__connection.commit()
        except Exception:
            logger.error("There was an error while selecting the transcript_id from the DB")
            raise Exception
        if record is not None:
            return record[0]
        else:
            return None


    def get_sample_id_by_name(self, sample_name, create_if_not_exists=False):
        # TODO: This function should cache results
        if create_if_not_exists:
            query = ("INSERT INTO samples(id, name) "
                     f"SELECT NULL,'{sample_name}' "
                     " WHERE NOT EXISTS ( "
                     "  SELECT 1"
                     "  FROM samples "
                     f"  WHERE name = '{sample_name}' "
                     ");"
                     )
            try:
                self.__cursor.execute(query)
            except Exception:
                logger.error("There was an error while inserting a new sample in the DB")
                raise Exception

        query = f"SELECT id from samples WHERE name = '{sample_name}'"
        try:
            self.__cursor.execute(query)
            record = self.__cursor.fetchone()
            self.__connection.commit()
        except Exception:
            logger.error("There was an error while selecting the sample_id from the DB")
            raise Exception
        if record is not None:
            return record[0]
        else:
            return None
