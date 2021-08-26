# -*- coding: utf-8 -*-

from enum import Enum
import datetime
import os
import sqlite3
import contextlib
from itertools import zip_longest, product

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

    table_defs = {} # table name -> column definitions (to be filled by derived classes)

    def __init__(self,
                 db_path:str,
                 create_mode=DBCreateMode.MUST_EXIST):
        self._db_path = db_path
        self._create_mode = create_mode
        self._connection = None
        self._cursor = None

    def _init_db(self):
        if self.table_defs:
            logger.debug("Setting up database tables")
            try:
                for table, column_defs in self.table_defs.items():
                    if type(column_defs) is not str: # list/tuple expected
                        column_defs = ", ".join(column_defs)
                    query = f"CREATE TABLE IF NOT EXISTS {table} ({column_defs})"
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

    # "reads" table:
    table_def_reads = ["id INTEGER NOT NULL PRIMARY KEY",
                       "name VARCHAR NOT NULL UNIQUE",
                       "sampleid INTEGER NOT NULL",
                       "transcriptid VARCHAR NOT NULL",
                       "refstart INT NOT NULL",
                       "refend INT NOT NULL",
                       "numevents INT NOT NULL",
                       "numsignals INT NOT NULL",
                       "dwelltime REAL NOT NULL",
                       "kmers INT NOT NULL",
                       "missing_kmers INT NOT NULL",
                       "NNNNN_kmers INT NOT NULL",
                       "mismatch_kmers INT NOT NULL",
                       "valid_kmers INT NOT NULL",
                       "FOREIGN KEY(sampleid) REFERENCES samples(id)",
                       "FOREIGN KEY(transcriptid) REFERENCES transcripts(id)"]

    # "kmer_sequences" table:
    table_def_kmer_seqs = ["id INTEGER NOT NULL PRIMARY KEY",
                           "sequence VARCHAR NOT NULL UNIQUE"]

    # "kmer_status" table:
    table_def_kmer_status = ["id INTEGER NOT NULL PRIMARY KEY",
                             "status VARCHAR NOT NULL UNIQUE"]

    # "kmers" table:
    # TODO: is combination of "readid" and "position" unique per kmer?
    # if so, use those as combined primary key (for access efficiency)?
    table_def_kmers = ["id INTEGER NOT NULL PRIMARY KEY",
                       "readid INTEGER NOT NULL",
                       "position INTEGER NOT NULL",
                       "sequenceid INTEGER",
                       # "sequence VARCHAR NOT NULL",
                       # "num_events INTEGER NOT NULL",
                       # "num_signals INTEGER NOT NULL",
                       "statusid INTEGER NOT NULL",
                       "dwell_time REAL NOT NULL",
                       # "NNNNN_dwell_time REAL NOT NULL",
                       # "mismatch_dwell_time REAL NOT NULL",
                       "median REAL NOT NULL",
                       "mad REAL NOT NULL",
                       "FOREIGN KEY(readid) REFERENCES reads(id)",
                       "FOREIGN KEY(sequenceid) REFERENCES kmer_sequences(id)",
                       "FOREIGN KEY(statusid) REFERENCES kmer_status(id)"]

    # "samples" table:
    table_def_samples = ["id INTEGER NOT NULL PRIMARY KEY",
                         "name VARCHAR NOT NULL UNIQUE",
                         "condition VARCHAR"]

    # "transcripts" table:
    table_def_transcripts = ["id INTEGER NOT NULL PRIMARY KEY",
                             "name VARCHAR NOT NULL UNIQUE"]

    table_defs = {"reads": table_def_reads,
                  "kmer_sequences": table_def_kmer_seqs,
                  "kmer_status": table_def_kmer_status,
                  "kmers": table_def_kmers,
                  "samples": table_def_samples,
                  "transcripts": table_def_transcripts}

    status_mapping = {"valid": 0, "NNNNN": 1, "mismatch": 2}
    sequence_mapping = {} # filled by "__init__"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        ## set up mapping table for sequences:
        self.sequence_mapping = {}
        seq_prod = product(["A", "C", "G", "T"], repeat=5)
        for i, seq in enumerate(seq_prod):
            self.sequence_mapping["".join(seq)] = i

    def _init_db(self):
        super()._init_db()
        ## fill "kmer_status" and "kmer_sequences" tables:
        self._cursor.executemany("INSERT INTO kmer_status VALUES (?, ?)",
                                 [(i, x) for x, i in self.status_mapping.items()])
        self._cursor.executemany("INSERT INTO kmer_sequences VALUES (?, ?)",
                                 [(i, x) for x, i in self.sequence_mapping.items()])
        self._connection.commit()

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
            status_id = self.status_mapping[res["status"]]
            # in case of unexpected kmer seq., this should give None (NULL in the DB):
            seq_id = self.sequence_mapping.get(res["ref_kmer"])
            self._cursor.execute("INSERT INTO kmers VALUES(NULL, ?, ?, ?, ?, ?, ?, ?)",
                                 (read_id, res["ref_pos"], seq_id, status_id,
                                  res["dwell_time"], res["median"], res["mad"]))
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

    # TODO: is this function never used?
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

    # "parameters" table:
    table_def_parameters = ["univariate_test VARCHAR CHECK (univariate_test in ('ST', 'MW', 'KS'))",
                            "gmm_covariance_type VARCHAR",
                            "gmm_test VARCHAR CHECK (gmm_test in ('anova', 'logit'))"]
    # TODO: add more parameters

    # "transcripts" table:
    table_def_transcripts = ["id INTEGER NOT NULL PRIMARY KEY",
                             "name VARCHAR NOT NULL UNIQUE"]

    # "whitelist" table:
    table_def_whitelist = ["transcriptid INTEGER NOT NULL",
                           "readid INTEGER NOT NULL UNIQUE", # foreign key for "reads" table in EventAlign DB
                           "FOREIGN KEY (transcriptid) REFERENCES transcripts(id)"]

    # "kmer_stats" table:
    table_def_kmer_stats = ["id INTEGER NOT NULL PRIMARY KEY",
                            "transcriptid INTEGER NOT NULL",
                            "kmer INTEGER NOT NULL",
                            "c1_mean_intensity REAL",
                            "c2_mean_intensity REAL",
                            "c1_median_intensity REAL",
                            "c2_median_intensity REAL",
                            "c1_sd_intensity REAL",
                            "c2_sd_intensity REAL",
                            "c1_mean_dwell REAL",
                            "c2_mean_dwell REAL",
                            "c1_median_dwell REAL",
                            "c2_median_dwell REAL",
                            "c1_sd_dwell REAL",
                            "c2_sd_dwell REAL",
                            "intensity_pvalue REAL",
                            "dwell_pvalue REAL",
                            "adj_intensity_pvalue REAL",
                            "adj_dwell_pvalue REAL",
                            "UNIQUE (transcriptid, kmer)",
                            "FOREIGN KEY (transcriptid) REFERENCES transcripts(id)"]
    # TODO: are "c1" and "c2" (conditions) properly defined?

    # "gmm_stats" table:
    table_def_gmm_stats = ["kmer_statsid INTEGER NOT NULL UNIQUE",
                           "n_components INTEGER NOT NULL",
                           "cluster_counts VARCHAR",
                           "test_stat REAL",
                           "test_pvalue REAL",
                           "adj_test_pvalue REAL",
                           "FOREIGN KEY (kmer_statsid) REFERENCES kmer_stats(id)"]

    table_defs = {"parameters": table_def_parameters,
                  "transcripts": table_def_transcripts,
                  "whitelist": table_def_whitelist,
                  "kmer_stats": table_def_kmer_stats}
    # table "gmm_stats" is only added when needed (see "__init__")

    def __init__(self,
                 db_path:str,
                 create_mode=DBCreateMode.MUST_EXIST,
                 with_gmm=True,
                 with_sequence_context=False):
        super().__init__(db_path, create_mode)
        self.__with_gmm = with_gmm
        self.__with_sequence_context = with_sequence_context
        if with_gmm:
            self.table_defs["gmm_stats"] = self.table_def_gmm_stats
        if with_sequence_context: # add additional columns for context p-values
            # column definitions must go before table constraints!
            constraints = self.table_defs["kmer_stats"][-2:]
            self.table_defs["kmer_stats"] = (self.table_defs["kmer_stats"][:-2] +
                                             ["intensity_pvalue_context REAL",
                                              "dwell_pvalue_context REAL",
                                              "adj_intensity_pvalue_context REAL",
                                              "adj_dwell_pvalue_context REAL"] +
                                             constraints)
            if with_gmm:
                constraints = self.table_defs["gmm_stats"][-1:]
                self.table_defs["gmm_stats"] = (self.table_defs["gmm_stats"][:-1] +
                                                ["test_pvalue_context REAL",
                                                 "adj_test_pvalue_context REAL"] +
                                                constraints)


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
        for kmer, res in test_results.items():
            values = [tx_id, kmer]
            values += res["shift_stats"].values()
            # insert 'None' (NULL) into adj. p-value columns:
            values += [res.get("intensity_pvalue"), res.get("dwell_pvalue"), None, None]
            if self.__with_sequence_context:
                values += [res.get("intensity_pvalue_context"), res.get("dwell_pvalue_context"), None, None]
            try:
                self._cursor.execute("INSERT INTO kmer_stats VALUES (NULL" + ", ?" * len(values) + ")", values)
            except:
                logger.error(f"Error storing statistics for transcript '{tx_name}', kmer {kmer}")
                raise
            kmer_statsid = self._cursor.lastrowid
            if self.__with_gmm:
                # insert 'None' (NULL) into adj. p-value columns:
                values = [kmer_statsid, res["gmm_model"].n_components, res.get("gmm_cluster_counts"),
                          res.get("gmm_test_stat"), res.get("gmm_pvalue"), None]
                if self.__with_sequence_context:
                    values += [res.get("gmm_pvalue_context"), None]
                qmarks = ", ".join(["?"] * len(values))
                try:
                    self._cursor.execute(f"INSERT INTO gmm_stats VALUES ({qmarks})", values)
                except:
                    logger.error(f"Error storing GMM stats for transcript '{tx_name}', kmer {kmer}")
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
