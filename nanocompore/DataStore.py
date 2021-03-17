# -*- coding: utf-8 -*-

from enum import Enum
import datetime
import os
import sqlite3
import contextlib

# Third party
from loguru import logger
from nanocompre.common import NanoporeError

class DataStore(object):
    """Store Nanocompore data in an SQLite database"""

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
                          "kmers INT NOT NULL,"
                          "missing_kmers INT NOT NULL,"
                          "NNNNN_kmers INT NOT NULL,"
                          "mismatch_kmers INT NOT NULL,"
                          "valid_kmers INT NOT NULL,"
                          "FOREIGN KEY(sampleid) REFERENCES samples(id)"
                          "FOREIGN KEY(transcriptid) REFERENCES transcripts(id)"
                          ")"
                          )

    create_kmers_query = ("CREATE TABLE IF NOT EXISTS kmers ("
                          "id INTEGER NOT NULL PRIMARY KEY,"
                          "readid INTEGER NOT NULL,"
                          "position INTEGER NOT NULL,"
                          "sequence INTEGER NOT NULL,"
                          "num_events INTEGER NOT NULL,"
                          "num_signals INTEGER NOT NULL,"
                          "status VARCHAR NOT NULL,"
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

    class DBCreateMode(Enum):
        """Options for handling (non-) existence of the SQLite database file"""
        MUST_EXIST = "r" # open for reading, error if file doesn't exist
        CREATE_MAYBE = "a" # use an existing database, otherwise create one
        OVERWRITE = "w" # always create a new database, overwrite if it exists


    def __init__(self,
                 db_path:str,
                 create_mode=DBCreateMode.MUST_EXIST):
        self.__db_path = db_path
        self.__create_mode = create_mode
        self.__connection = None
        self.__cursor = None

    def __enter__(self):
        if self.__create_mode == DBCreateMode.MUST_EXIST and not os.path.exists(self.__db_path):
            raise NanocomporeError(f"Database file '{self.__db_path}' does not exist")
        if self.__create_mode == DBCreateMode.OVERWRITE:
            with contextlib.suppress(FileNotFoundError): # file may not exist
                os.remove(self.__db_path)
                logger.debug(f"Removed existing database file '{self.__db_path}'")
        try:
            logger.debug("Connecting to database")
            self.__connection = sqlite3.connect(self.__db_path)
            self.__connection.row_factory = sqlite3.Row
            self.__cursor = self.__connection.cursor()
        except:
            logger.error("Error connecting to database")
            raise
        if self.__create_mode == DBCreateMode.OVERWRITE or \
           (self.__create_mode == DBCreateMode.CREATE_MAYBE and not os.path.exists(self.__db_path)):
            self.__init_db()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.__connection:
            logger.debug("Closing database connection")
            self.__connection.commit()
            self.__connection.close()
            self.__connection = None
            self.__cursor = None

   def __init_db(self):
        logger.debug("Setting up database tables")
        try:
            self.__cursor.execute(self.create_reads_query)
            self.__cursor.execute(self.create_kmers_query)
            self.__cursor.execute(self.create_samples_query)
            self.__cursor.execute(self.create_transcripts_query)
            self.__connection.commit()
        except:
            logger.error("Error creating database tables")
            raise

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
        values = (read.read_id, sample_id, tx_id, read.ref_start, read.ref_end,
                  read.n_events, read.n_signals, read.dwell_time) + tuple(read.kmers_status.values())
        try:
            self.__cursor.execute("INSERT INTO reads VALUES(NULL" + ", ?" * len(values) + ")",
                                  values)
            read_id = self.__cursor.lastrowid
        except Exception:
            logger.error("Error inserting read into database")
            raise Exception

        for kmer in read.kmer_l:
            self.__store_kmer(kmer=kmer, read_id=read_id)
        self.__connection.commit()
        # TODO check for success and return true/false

    def __store_kmer(self, kmer, read_id):
        """
        Insert data corresponding to a kmer into the DB.
        Args:
            kmer (Kmer): instance of class Kmer that contains kmer-level data
            read_id (int): DB key of the read the kmer belongs to
        """
        res = kmer.get_results() # needed for 'median' and 'mad' values
        try:
            self.__cursor.execute("INSERT INTO kmers VALUES(NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                              (read_id, res["ref_pos"], res["ref_kmer"], res["num_events"],
                               res["num_signals"], res["status"], res["dwell_time"],
                               res["NNNNN_dwell_time"], res["mismatch_dwell_time"], res["median"], res["mad"]))
        except Exception:
            logger.error("Error inserting kmer into database")
            raise Exception

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
                logger.error("Error while inserting transcript into the database")
                raise Exception

        query = f"SELECT id from transcripts WHERE name = '{tx_name}'"
        try:
            self.__cursor.execute(query)
            record = self.__cursor.fetchone()
            self.__connection.commit()
        except Exception:
            logger.error("Error while selecting transcript ID from the database")
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
                logger.error("Error while inserting sample into the database")
                raise Exception

        query = f"SELECT id from samples WHERE name = '{sample_name}'"
        try:
            self.__cursor.execute(query)
            record = self.__cursor.fetchone()
            self.__connection.commit()
        except Exception:
            logger.error("Error while selecting sample ID from the database")
            raise Exception
        if record is not None:
            return record[0]
        else:
            return None

    @property
    def cursor(self):
        return self.__cursor

    def get_samples(self, sample_dict=None):
        if not self.__connection:
            raise NanocomporeError("Database connection not yet opened")
        expected_samples = []
        if sample_dict: # query only relevant samples
            for samples in sample_dict.values():
                expected_samples += samples
            if not expected_samples:
                raise NanocomporeError("No sample names in 'sample_dict'")
            where = " WHERE name IN ('%s')" % "', '".join(expected_samples)
        else:
            where = ""
        db_samples = {}
        try:
            self.__cursor.execute("SELECT * FROM samples" + where)
            for row in self.__cursor:
                db_samples[row["id"]] = row["name"]
        except Exception:
            logger.error("Error reading sample names from database")
            raise Exception
        for sample in expected_samples: # check that requested samples are in DB
            if sample not in db_samples.values():
                raise NanocomporeError(f"Sample '{sample}' not present in database")
        return db_samples
