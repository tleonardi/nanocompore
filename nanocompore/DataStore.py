# -*- coding: utf-8 -*-

from enum import Enum
import datetime
import os
import sqlite3
import contextlib
from itertools import zip_longest

# Third party
from loguru import logger
from nanocompore.common import NanocomporeError


class DBCreateMode(Enum):
    """Options for handling (non-) existence of the SQLite database file"""
    MUST_EXIST = "r" # open for reading, error if file doesn't exist
    CREATE_MAYBE = "a" # use an existing database, otherwise create one
    OVERWRITE = "w" # always create a new database, overwrite if it exists


class DataStore(object):
    """Store Nanocompore data in an SQLite database - base class"""

    create_tables_queries = {} # table name -> SQL query (to be filled by derived classes)

    def __init__(self,
                 db_path:str,
                 create_mode=DBCreateMode.MUST_EXIST):
        self._db_path = db_path
        self._create_mode = create_mode
        self._connection = None
        self._cursor = None

    def _init_db(self):
        if self.create_tables_queries:
            logger.debug("Setting up database tables")
            try:
                for table, query in self.create_tables_queries.items():
                    self._cursor.execute(query)
                self._connection.commit()
            except:
                logger.error(f"Error creating database table '{table}'")
                raise

    def __enter__(self):
        init_db = False
        if self._create_mode == DBCreateMode.MUST_EXIST and not os.path.exists(self._db_path):
            raise NanocomporeError(f"Database file '{self._db_path}' does not exist")
        if self._create_mode == DBCreateMode.OVERWRITE:
            with contextlib.suppress(FileNotFoundError): # file may not exist
                os.remove(self._db_path)
                logger.debug(f"Removed existing database file '{self._db_path}'")
            init_db = True
        if self._create_mode == DBCreateMode.CREATE_MAYBE and not os.path.exists(self._db_path):
            init_db = True
        try:
            logger.debug("Connecting to database")
            self._connection = sqlite3.connect(self._db_path)
            self._connection.row_factory = sqlite3.Row
            self._cursor = self._connection.cursor()
        except:
            logger.error("Error connecting to database")
            raise
        if init_db:
            self._init_db()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._connection:
            logger.debug("Closing database connection")
            self._connection.commit()
            self._connection.close()
            self._connection = None
            self._cursor = None

    @property
    def cursor(self):
        return self._cursor


class DataStore_EventAlign(DataStore):
    """Store Nanocompore data in an SQLite database - subclass for Eventalign_collapse results"""

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
    # TODO: encode 'status' as int to save space (foreign key referencing a table with all possible statuses)

    create_samples_query = ("CREATE TABLE IF NOT EXISTS samples ("
                            "id INTEGER NOT NULL PRIMARY KEY,"
                            "name VARCHAR NOT NULL UNIQUE,"
                            "condition VARCHAR"
                            ")"
                            )

    create_transcripts_query = ("CREATE TABLE IF NOT EXISTS transcripts ("
                                "id INTEGER NOT NULL PRIMARY KEY,"
                                "name VARCHAR NOT NULL UNIQUE"
                                ")"
                                )

    create_tables_queries = {"reads": create_reads_query,
                             "kmers": create_kmers_query,
                             "samples": create_samples_query,
                             "transcripts": create_transcripts_query}

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
            self._cursor.execute("INSERT INTO reads VALUES(NULL" + ", ?" * len(values) + ")",
                                  values)
            read_id = self._cursor.lastrowid
        except Exception:
            logger.error("Error inserting read into database")
            raise Exception

        for kmer in read.kmer_l:
            self.__store_kmer(kmer=kmer, read_id=read_id)
        self._connection.commit()
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
            self._cursor.execute("INSERT INTO kmers VALUES(NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
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
                self._cursor.execute(query)
            except Exception:
                logger.error("Error while inserting transcript into the database")
                raise Exception

        query = f"SELECT id from transcripts WHERE name = '{tx_name}'"
        try:
            self._cursor.execute(query)
            record = self._cursor.fetchone()
            self._connection.commit()
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
                self._cursor.execute(query)
            except Exception:
                logger.error("Error while inserting sample into the database")
                raise Exception

        query = f"SELECT id from samples WHERE name = '{sample_name}'"
        try:
            self._cursor.execute(query)
            record = self._cursor.fetchone()
            self._connection.commit()
        except Exception:
            logger.error("Error while selecting sample ID from the database")
            raise Exception
        if record is not None:
            return record[0]
        else:
            return None

    def get_samples(self, sample_dict=None):
        if not self._connection:
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
            self._cursor.execute("SELECT * FROM samples" + where)
            for row in self._cursor:
                db_samples[row["id"]] = row["name"]
        except Exception:
            logger.error("Error reading sample names from database")
            raise Exception
        for sample in expected_samples: # check that requested samples are in DB
            if sample not in db_samples.values():
                raise NanocomporeError(f"Sample '{sample}' not present in database")
        return db_samples

    def store_sample_info(self, sample_dict):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        # query: insert sample; if it exists, update condition if that's missing
        query = "INSERT INTO samples(id, name, condition) VALUES (NULL, ?, ?) " \
            "ON CONFLICT(name) DO UPDATE SET condition = excluded.condition " \
            "WHERE condition IS NULL"
        for condition, samples in sample_dict.items():
            try:
                self._cursor.executemany(query, [(condition, sample) for sample in samples])
            except:
                logger.error(f"Error storing sample information for condition '{condition}'")
                raise
        self._connection.commit()


class DataStore_SampComp(DataStore):
    """Store Nanocompore data in an SQLite database - subclass for SampComp results"""

    create_transcripts_query = ("CREATE TABLE IF NOT EXISTS transcripts ("
                                "id INTEGER NOT NULL PRIMARY KEY,"
                                "name VARCHAR NOT NULL UNIQUE"
                                ")")

    create_whitelist_query = ("CREATE TABLE IF NOT EXISTS whitelist ("
                              "transcriptid INTEGER NOT NULL,"
                              "readid INTEGER NOT NULL UNIQUE,"
                              "FOREIGN KEY (transcriptid) REFERENCES transcripts(id)"
                              # "readid" is foreign key for "reads" table in EventAlign DB
                              ")")

    create_kmer_stats_query = ("CREATE TABLE IF NOT EXISTS kmer_stats ("
                               "id INTEGER NOT NULL PRIMARY KEY,"
                               "transcriptid INTEGER NOT NULL,"
                               "kmer INTEGER NOT NULL,"
                               "c1_mean_intensity REAL,"
                               "c2_mean_intensity REAL,"
                               "c1_median_intensity REAL,"
                               "c2_median_intensity REAL,"
                               "c1_sd_intensity REAL,"
                               "c2_sd_intensity REAL,"
                               "c1_mean_dwell REAL,"
                               "c2_mean_dwell REAL,"
                               "c1_median_dwell REAL,"
                               "c2_median_dwell REAL,"
                               "c1_sd_dwell REAL,"
                               "c2_sd_dwell REAL,"
                               "UNIQUE (transcriptid, kmer),"
                               "FOREIGN KEY (transcriptid) REFERENCES transcripts(id)"
                               ")")
    # TODO: are "c1" and "c2" (conditions) properly defined?

    create_gmm_stats_query = ("CREATE TABLE IF NOT EXISTS gmm_stats ("
                              "kmer_statsid INTEGER NOT NULL UNIQUE,"
                              "covariance_type VARCHAR,"
                              "n_components INTEGER,"
                              "cluster_counts VARCHAR,"
                              "FOREIGN KEY (kmer_statsid) REFERENCES kmer_stats(id)"
                              ")")
    # TODO: store GMM cluster counts in a separate table (one row per sample)
    # TODO: if "covariance_type" is the same for all rows, store it in a "parameters" table?

    # TODO: add column for adjusted p-values in tables below?
    create_gmm_results_query = ("CREATE TABLE IF NOT EXISTS gmm_results ("
                                "gmm_statsid INTEGER NOT NULL UNIQUE,"
                                "test VARCHAR NOT NULL CHECK (test in ('anova', 'logit')),"
                                "test_pvalue REAL,"
                                "test_stat REAL," # anova: delta logit, logit: LOR
                                "UNIQUE (gmm_statsid, test),"
                                "FOREIGN KEY (gmm_statsid) REFERENCES gmm_stats(id)"
                                ")")

    create_univariate_results_query = ("CREATE TABLE IF NOT EXISTS univariate_results ("
                                       "kmer_statsid INTEGER NOT NULL,"
                                       "test VARCHAR NOT NULL CHECK (test in ('ST', 'MW', 'KS')),"
                                       "intensity_pvalue REAL,"
                                       "dwell_pvalue REAL,"
                                       "UNIQUE (kmer_statsid, test),"
                                       "FOREIGN KEY (kmer_statsid) REFERENCES kmer_stats(id)"
                                       ")")

    create_tables_queries = {"transcripts": create_transcripts_query,
                             "whitelist": create_whitelist_query,
                             "kmer_stats": create_kmer_stats_query,
                             "gmm_stats": create_gmm_stats_query,
                             "gmm_results": create_gmm_results_query,
                             "univariate_results": create_univariate_results_query}

    def __insert_transcript_get_id(self, tx_name):
        try:
            self._cursor.execute("SELECT id FROM transcripts WHERE name = ?", [tx_name])
            if (row := self._cursor.fetchone()) is not None:
                return row["id"]
            self._cursor.execute("INSERT INTO transcripts VALUES (NULL, ?)", [tx_name])
            self._connection.commit()
            # TODO: if there could be multiple writing threads, "INSERT OR IGNORE"
            # query should go before "SELECT"
            return self._cursor.lastrowid
        except:
            logger.error(f"Failed to insert/look up transcript '{tx_name}'")
            raise


    def store_test_results(self, tx_name, test_results):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        tx_id = self.__insert_transcript_get_id(tx_name)
        univar_pvalues = [f"{t}_{m}_pvalue" for t in ["MW", "KS", "ST"]
                          for m in ["intensity", "dwell"]]
        for kmer, res in test_results.items():
            values = [tx_id, kmer]
            values += res["shift_stats"].values()
            try:
                self._cursor.execute("INSERT INTO kmer_stats VALUES (NULL" + ", ?" * len(values) + ")", values)
            except:
                logger.error(f"Error storing statistics for transcript '{tx_name}', kmer {kmer}")
                raise
            kmer_statsid = self._cursor.lastrowid
            for test in ["MW", "KS", "ST"]:
                ipv = res.get(test + "_intensity_pvalue")
                dpv = res.get(test + "_dwell_pvalue")
                if (ipv is not None) or (dpv is not None): # can't use ':=' here because we need both values
                    try:
                        self._cursor.execute("INSERT INTO univariate_results VALUES (?, ?, ?, ?)",
                                             (kmer_statsid, test, ipv, dpv))
                    except:
                        logger.error(f"Error storing {test} test results for transcript '{tx_name}', kmer {kmer}")
                        raise
            if "GMM_model" in res:
                try:
                    self._cursor.execute("INSERT INTO gmm_stats VALUES (?, ?, ?, ?)",
                                         (kmer_statsid,
                                          res["GMM_model"]["model"].covariance_type,
                                          res["GMM_model"]["model"].n_components,
                                          res["GMM_model"]["cluster_counts"]))
                except:
                    logger.error(f"Error storing GMM stats for transcript '{tx_name}', kmer {kmer}")
                    raise
                gmm_statsid = self._cursor.lastrowid
                # store results of logit and/or ANOVA test on GMM components:
                test_stats = {"logit": "coef", "anova": "delta_logit"}
                for test, stat in test_stats.items():
                    if f"GMM_{test}_pvalue" in res:
                        try:
                            self._cursor.execute("INSERT INTO gmm_results VALUES (?, ?, ?, ?)",
                                                 (gmm_statsid, test,
                                                  res[f"GMM_{test}_pvalue"],
                                                  res[f"GMM_{test}_model"][stat]))
                        except:
                            logger.error(f"Error storing GMM {test} results for transcript '{tx_name}', kmer {kmer}")
                            raise
            self._connection.commit()


    def store_whitelist(self, whitelist):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        for tx_name, read_dict in whitelist:
            try:
                tx_id = self.__insert_transcript_get_id(tx_name)
                for cond_reads in read_dict.values():
                    for sample_reads in cond_reads.values():
                        values = zip_longest([], sample_reads, fillvalue=tx_id)
                        self._cursor.executemany("INSERT INTO whitelist VALUES (?, ?)", values)
                        # TODO: store sample/condition information (again)?
                        # it can be retrieved from "reads"/"samples" tables given "readid"
                self._connection.commit()
            except:
                logger.error(f"Error storing whitelisted reads for transcript '{tx_name}'")
                raise
