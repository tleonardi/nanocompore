# -*- coding: utf-8 -*-

from enum import Enum
import datetime
import os
import sqlite3
import contextlib
import math
from itertools import zip_longest, product

# Third party
import numpy as np
from loguru import logger
from sklearn.mixture import GaussianMixture
from nanocompore.common import NanocomporeError


class DBCreateMode(Enum):
    """Options for handling (non-) existence of the SQLite database file"""
    MUST_EXIST = "r" # open for reading, error if file doesn't exist
    CREATE_MAYBE = "a" # use an existing database, otherwise create one
    OVERWRITE = "w" # always create a new database, overwrite if it exists


class DataStore:
    """Store Nanocompore data in an SQLite database - base class"""

    table_defs = {} # table name -> column definitions (to be filled by derived classes)
    status_mapping = {"valid": 0, "NNNNN": 1, "mismatch": 2}
    sequence_mapping = {} # filled by "__init__"

    def __init__(self,
                 db_path:str,
                 create_mode=DBCreateMode.MUST_EXIST):
        self._db_path = db_path
        self._create_mode = create_mode
        self._connection = None
        self._cursor = None
        self.sequence_mapping = {}
        seq_prod = product(["A", "C", "G", "T"], repeat=5)
        for i, seq in enumerate(seq_prod):
            self.sequence_mapping["".join(seq)] = i

    def _create_tables(self, table_defs):
        try:
            for table, column_defs in table_defs.items():
                if type(column_defs) is not str: # list/tuple expected
                    column_defs = ", ".join(column_defs)
                query = f"CREATE TABLE IF NOT EXISTS {table} ({column_defs})"
                self._cursor.execute(query)
            self._connection.commit()
        except:
            logger.error(f"Error creating database table '{table}'")
            raise

    def _init_db(self):
        if self.table_defs:
            logger.debug("Setting up database tables")
            self._create_tables(self.table_defs)

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
            logger.trace("Connecting to database")
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
            logger.trace("Closing database connection")
            self._connection.commit()
            self._connection.close()
            self._connection = None
            self._cursor = None

    @property
    def cursor(self):
        return self._cursor

    def add_or_reset_column(self, table, col_name, col_def, reset_value=None):
        """Try to add a column to a table. If the column already exists, reset its values."""
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        try:
            self._cursor.execute(f"ALTER TABLE {table} ADD COLUMN {col_name} {col_def}")
        except sqlite3.OperationalError as error:
            if error.args[0].startswith("duplicate column name:"):
                self._cursor.execute(f"UPDATE {table} SET {col_name} = ?", (reset_value, ))
            else:
                logger.error("Error adding column '{col_name}' to table '{table}'")
                raise
        except:
            logger.error("Error adding column '{col_name}' to table '{table}'")
            raise
        self._connection.commit()


class DataStore_master(DataStore):
    """Store Nanocompore data in an SQLite database - master database"""

    # TODO: add "parameters" table and store parameters

    # "kmer_sequences" table:
    table_def_kmer_seqs = ["id INTEGER NOT NULL PRIMARY KEY",
                           "sequence VARCHAR NOT NULL UNIQUE"]

    # "kmer_status" table:
    table_def_kmer_status = ["id INTEGER NOT NULL PRIMARY KEY",
                             "status VARCHAR NOT NULL UNIQUE"]

    # "transcript_status" table:
    table_def_transcript_status = ["id INTEGER NOT NULL PRIMARY KEY",
                                   "status VARCHAR NOT NULL UNIQUE"]

    # "samples" table:
    table_def_samples = ["id INTEGER NOT NULL PRIMARY KEY",
                         "name VARCHAR NOT NULL UNIQUE",
                         "file VARCHAR NOT NULL UNIQUE",
                         "condition VARCHAR NOT NULL CHECK (condition != '')"]

    # "transcripts" table:
    table_def_transcripts = ["id INTEGER NOT NULL PRIMARY KEY",
                             "name VARCHAR NOT NULL UNIQUE",
                             "subdir VARCHAR NOT NULL",
                             "status INTEGER",
                             "FOREIGN KEY(status) REFERENCES transcript_status(id)"]

    # "parameters" table:
    table_def_parameters = ["step VARCHAR NOT NULL",
                            "name VARCHAR NOT NULL",
                            "value VARCHAR",
                            "type VARCHAR NOT NULL",
                            "UNIQUE (step, name)"]

    table_defs = {"kmer_sequences": table_def_kmer_seqs,
                  "kmer_status": table_def_kmer_status,
                  "transcript_status": table_def_transcript_status,
                  "samples": table_def_samples,
                  "transcripts": table_def_transcripts,
                  "parameters": table_def_parameters}

    def _init_db(self):
        super()._init_db()
        ## fill "kmer_status", "kmer_sequences" and "transcript_status" tables:
        self._cursor.executemany("INSERT INTO kmer_status VALUES (?, ?)",
                                 [(i, x) for x, i in self.status_mapping.items()])
        self._cursor.executemany("INSERT INTO kmer_sequences VALUES (?, ?)",
                                 [(i, x) for x, i in self.sequence_mapping.items()])
        tx_statuses = [(0, "valid"), (1, "too short"), (2, "low coverage"), (3, "no data")]
        self._cursor.executemany("INSERT INTO transcript_status VALUES (?, ?)", tx_statuses)
        self._connection.commit()

    def store_sample(self, sample_name, file_path, condition):
        """Store a new sample in the database"""
        try:
            self._cursor.execute("INSERT INTO samples VALUES(NULL, ?, ?, ?)",
                                 (sample_name, file_path, condition))
            self._connection.commit()
            return self._cursor.lastrowid
        except Exception:
            logger.error("Error inserting sample into database")
            raise Exception

    def store_transcript(self, name, subdir=""):
        """Store a transcript in the database (if not already stored)"""
        try:
            self._cursor.execute("SELECT id FROM transcripts WHERE name = ?", [name])
            if (row := self._cursor.fetchone()) is not None:
                return row["id"]
            self._cursor.execute("INSERT INTO transcripts(name, subdir) VALUES (?, ?)", (name, subdir))
            self._connection.commit()
            # TODO: if there could be multiple writing threads, "INSERT OR IGNORE"
            # query should go before "SELECT"
            return self._cursor.lastrowid
        except:
            logger.error(f"Failed to insert transcript '{name}'")
            raise

    def get_sample_info(self):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        db_samples = {}
        try:
            self._cursor.execute("SELECT * FROM samples ORDER BY id")
            for row in self._cursor:
                # TODO: include file in output?
                db_samples.setdefault(row["condition"], []).append((row["id"], row["name"]))
        except:
            logger.error("Error reading sample information from database")
            raise
        return db_samples

    def init_test_results(self, univariate_test=True, gmm_test=True,
                          sequence_context=False, drop_old=True):
        if not univariate_test and not gmm_test: # no test results to store
            return
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        # store settings that impact database schema:
        self._univariate_test = univariate_test
        self._gmm_test = gmm_test
        self._sequence_context = sequence_context
        # TODO: add "id" column?
        table_def = ["transcriptid INTEGER NOT NULL",
                     "kmer_pos INTEGER NOT NULL"]
        # add more table columns as necessary:
        if univariate_test:
            self.add_or_reset_column("transcripts", "n_univariate_tests", "INTEGER")
            table_def += ["intensity_pvalue REAL", "dwelltime_pvalue REAL"]
            if sequence_context:
                table_def += ["intensity_pvalue_context REAL", "dwelltime_pvalue_context REAL"]
        if gmm_test:
            self.add_or_reset_column("transcripts", "n_gmm_tests", "INTEGER")
            table_def.append("gmm_pvalue REAL")
            if sequence_context:
                table_def.append("gmm_pvalue_context REAL")
        table_def.append("FOREIGN KEY(transcriptid) REFERENCES transcripts(id)")
        table_def = ", ".join(table_def)
        try:
            if drop_old: # remove previous results if any exist
                self._cursor.execute("DROP TABLE IF EXISTS test_results")
            self._cursor.execute(f"CREATE TABLE IF NOT EXISTS test_results ({table_def})")
        except:
            logger.error("Error creating 'test_results' table")
            raise
        self._connection.commit()

    def store_test_results(self, tx_id, status, results):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        # store transcript status and numbers of tests performed (for multiple testing correction):
        assign = [f"status = {status}"]
        if self._univariate_test and results:
            assign.append("n_univariate_tests = %d" % results["n_univariate_tests"])
        if self._gmm_test and results:
            assign.append("n_gmm_tests = %d" % results["n_gmm_tests"])
        assign = ", ".join(assign)
        sql = f"UPDATE transcripts SET {assign} WHERE id = {tx_id}"
        try:
            self._cursor.execute(sql)
            self._connection.commit()
        except:
            logger.error(f"Error updating status and test counts for transcript {tx_id}")
            raise
        # store kmers with significant test results:
        # TODO: include kmer sequence?
        if results and results["test_results"]:
            for kmer, res in results["test_results"].items():
                values = [tx_id, kmer]
                if self._univariate_test:
                    values += [res["intensity_pvalue"], res["dwelltime_pvalue"]]
                    if self._sequence_context:
                        values += [res["intensity_pvalue_context"], res["dwelltime_pvalue_context"]]
                if self._gmm_test:
                    values.append(res["gmm_pvalue"])
                    if self._sequence_context:
                        values.append(res["gmm_pvalue_context"])
                qmarks = ", ".join(["?"] * len(values))
                try:
                    self._cursor.execute(f"INSERT INTO test_results VALUES ({qmarks})", values)
                except:
                    logger.error(f"Error storing statistics for kmer {kmer}")
                    raise
        self._connection.commit()

    def store_parameters(self, step, **kwargs):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        values = [(step, k, str(v), type(v).__name__) for k, v in kwargs.items()]
        sql = "INSERT INTO parameters(step, name, value, type) VALUES (?, ?, ?, ?) ON CONFLICT(step, name) DO UPDATE SET value = excluded.value, type = excluded.type"
        try:
            self._cursor.executemany(sql, values)
            self._connection.commit()
        except:
            logger.error(f"Error storing parameters")
            raise

    def reset_SampComp(self):
        """Reset the database to a state without SampComp results"""
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        self._cursor.execute("DELETE FROM parameters WHERE step = 'SC'")
        self._cursor.execute("DROP TABLE IF EXISTS test_results")
        # "ALTER table DROP COLUMN" not supported by SQLite version used here
        # (sqlite3.sqlite_version: '3.31.1' - supported from 3.35)
        # can't remove columns, so reset values:
        for column in ("n_univariate_tests", "n_gmm_tests"):
            try:
                self._cursor.execute(f"UPDATE transcripts SET {column} = NULL")
            except sqlite3.OperationalError as error:
                # ignore "no such column" error:
                if not error.args[0].startswith("no such column:"):
                    raise
        self._connection.commit()


class DataStore_transcript(DataStore):
    """Store Nanocompore data in an SQLite database - single-transcript database"""

    # "transcript" table (information is also stored in master DB):
    table_def_transcript = ["id INTEGER NOT NULL PRIMARY KEY",
                            "name VARCHAR NOT NULL UNIQUE"]
    # TODO: store transcript sequence here instead of kmer seqs. in "kmers" table?

    # "reads" table:
    table_def_reads = ["id INTEGER NOT NULL PRIMARY KEY",
                       "name VARCHAR NOT NULL UNIQUE",
                       "sampleid INTEGER NOT NULL", # references 'samples(id)' in master DB
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
                       "intensity_mean REAL NOT NULL",
                       "intensity_sd REAL NOT NULL",
                       "dwelltime_log10_mean REAL NOT NULL",
                       "dwelltime_log10_sd REAL NOT NULL"]

    # "kmers" table:
    table_def_kmers = ["readid INTEGER NOT NULL",
                       "position INTEGER NOT NULL",
                       "sequenceid INTEGER", # references 'kmer_sequences(id)' in master DB
                       # "num_events INTEGER NOT NULL",
                       # "num_signals INTEGER NOT NULL",
                       "statusid INTEGER NOT NULL", # references 'kmer_status(id)' in master DB
                       "dwelltime_log10 REAL NOT NULL",
                       # "NNNNN_dwelltime_log10 REAL NOT NULL",
                       # "mismatch_dwelltime_log10 REAL NOT NULL",
                       "intensity REAL NOT NULL",
                       # "intensity_deviation REAL NOT NULL",
                       "PRIMARY KEY (readid, position)",
                       "FOREIGN KEY (readid) REFERENCES reads(id)"]

    # "kmer_stats" table:
    table_def_kmer_stats = ["kmer_pos INTEGER NOT NULL PRIMARY KEY",
                            # stats are after scaling, i.e. based on z-scores:
                            "c1_mean_intensity REAL",
                            "c2_mean_intensity REAL",
                            "c1_median_intensity REAL",
                            "c2_median_intensity REAL",
                            "c1_sd_intensity REAL",
                            "c2_sd_intensity REAL",
                            "c1_mean_dwelltime REAL",
                            "c2_mean_dwelltime REAL",
                            "c1_median_dwelltime REAL",
                            "c2_median_dwelltime REAL",
                            "c1_sd_dwelltime REAL",
                            "c2_sd_dwelltime REAL",
                            "intensity_pvalue REAL",
                            "dwelltime_pvalue REAL"]
    # TODO: are "c1" and "c2" (conditions) properly defined?

    # "gmm_stats" table:
    table_def_gmm_stats = ["kmer_pos INTEGER NOT NULL PRIMARY KEY",
                           "cluster_counts VARCHAR",
                           "test_stat REAL",
                           "test_pvalue REAL",
                           "FOREIGN KEY (kmer_pos) REFERENCES kmer_stats(kmer_pos)"]

    # "gmm_components" table:
    # TODO: store other 'GaussianMixture' attributes (e.g. 'converged_', 'n_iter_')?
    # TODO: add 'coveriance_type' if it becomes variable (currently fixed to 'full' in TxComp)
    # TODO: split into separate 'gmm_models' and 'gmm_components' tables to avoid redundancy?
    table_def_gmm_components = ["kmer_pos INTEGER NOT NULL",
                                "component INTEGER NOT NULL",
                                "weight REAL",
                                # stats are after scaling, i.e. based on z-scores:
                                "mean_intensity REAL NOT NULL",
                                "mean_dwelltime REAL NOT NULL",
                                # elements of the Cholesky decomposition of the
                                # precision matrix (inverse of covariance matrix);
                                # matrix is triangular, so (1,0) is always 0:
                                "precision_cholesky_00 REAL NOT NULL",
                                "precision_cholesky_01 REAL NOT NULL",
                                "precision_cholesky_11 REAL NOT NULL",
                                "model_bic REAL",
                                "PRIMARY KEY (kmer_pos, component)",
                                "FOREIGN KEY (kmer_pos) REFERENCES kmer_stats(kmer_pos)"]

    table_defs = {"reads": table_def_reads,
                  "kmers": table_def_kmers,
                  "transcript": table_def_transcript}
    # "..._stats" and "gmm_components" tables are added later, if needed

    def __init__(self, db_path:str, tx_name:str, tx_id:int, create_mode=DBCreateMode.MUST_EXIST):
        self.tx_name = tx_name
        self.tx_id = tx_id
        super().__init__(db_path, create_mode)

    def _init_db(self):
        super()._init_db()
        self._cursor.execute("INSERT INTO transcript VALUES (?, ?)", (self.tx_id, self.tx_name))
        self._connection.commit()

    def store_read(self, read, sample_id):
        """
        Insert data corresponding to a new read into the DB.
        Args:
            read (Read): an instance of class Read that contains read-level data, incl. kmer data.
            sample_id: DB ID of the sample to which the read belongs
        """
        values = (read.read_id, sample_id, read.ref_start, read.ref_end,
                  read.n_events, read.n_signals, read.dwell_time) + \
                  tuple(read.kmers_status.values()) + \
                  tuple(read.get_kmer_stats().values())
        try:
            self._cursor.execute("INSERT INTO reads VALUES(NULL" + ", ?" * len(values) + ")",
                                  values)
            read_id = self._cursor.lastrowid
        except Exception:
            logger.error("Error inserting read into database")
            raise Exception

        for kmer in read.kmer_l:
            self._store_kmer(kmer, read_id)
        self._connection.commit()

    def _store_kmer(self, kmer, read_id):
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
            self._cursor.execute("INSERT INTO kmers VALUES(?, ?, ?, ?, ?, ?)",
                                 (read_id, res["ref_pos"], seq_id, status_id,
                                  math.log10(res["dwell_time"]), res["median"]))
        except Exception:
            logger.error("Error inserting kmer into database")
            raise Exception

    def create_stats_tables(self, with_gmm=True, with_gmm_components="none",
                            with_sequence_context=False, drop_old=True):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        self._with_gmm = with_gmm
        self._with_gmm_components = with_gmm_components
        self._with_sequence_context = with_sequence_context
        table_defs = {"kmer_stats": self.table_def_kmer_stats}
        if with_gmm:
            table_defs["gmm_stats"] = self.table_def_gmm_stats
            if with_gmm_components != "none":
                table_defs["gmm_components"] = self.table_def_gmm_components
        if with_sequence_context: # add additional columns for context p-values
            table_defs["kmer_stats"] += ["intensity_pvalue_context REAL",
                                         "dwelltime_pvalue_context REAL"]
            if with_gmm:
                # column definitions must go before table constraints!
                constraint = table_defs["gmm_stats"].pop()
                table_defs["gmm_stats"] += ["test_pvalue_context REAL", constraint]
        if drop_old: # remove previous results if any exist
            self._cursor.execute("DROP TABLE IF EXISTS gmm_components")
            self._cursor.execute("DROP TABLE IF EXISTS gmm_stats")
            self._cursor.execute("DROP TABLE IF EXISTS kmer_stats")
        self._create_tables(table_defs)

    def store_test_results(self, test_results):
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        for kmer, res in test_results.items():
            values = [kmer]
            values += res["shift_stats"].values()
            values += [res.get("intensity_pvalue"), res.get("dwelltime_pvalue")]
            if self._with_sequence_context:
                values += [res.get("intensity_pvalue_context"), res.get("dwelltime_pvalue_context")]
            qmarks = ", ".join(["?"] * len(values))
            try:
                self._cursor.execute(f"INSERT INTO kmer_stats VALUES ({qmarks})", values)
            except:
                logger.error(f"Error storing statistics for kmer {kmer}")
                raise
            if self._with_gmm:
                best_index = res["gmm_best_index"]
                if self._with_gmm_components == "best": # store components of best GMM
                    self.store_GMM_components(res["gmm_models"][best_index], res["gmm_bics"][best_index], kmer)
                elif self._with_gmm_components == "all": # store components of all GMMs
                    for index in range(len(res["gmm_models"])):
                        self.store_GMM_components(res["gmm_models"][index], res["gmm_bics"][index], kmer)
                # for GMM stats, skip uninformative GMMs with only one component:
                if res["gmm_models"][best_index].n_components > 1:
                    values = [kmer, res.get("gmm_cluster_counts"), res.get("gmm_test_stat"), res.get("gmm_pvalue")]
                    if self._with_sequence_context:
                        values += [res.get("gmm_pvalue_context")]
                    qmarks = ", ".join(["?"] * len(values))
                    try:
                        self._cursor.execute(f"INSERT INTO gmm_stats VALUES ({qmarks})", values)
                    except:
                        logger.error(f"Error storing GMM stats for kmer {kmer}")
                        raise
            self._connection.commit()

    def store_GMM_components(self, gmm, bic, kmer_pos):
        for i in range(gmm.n_components):
            # use index 0 for 1-component GMM, indexes 1/2 for 2-component GMMs:
            index = i + int(gmm.n_components > 1)
            values = [kmer_pos, index, gmm.weights_[i]] + gmm.means_[i].tolist()
            cholesky = gmm.precisions_cholesky_[i].flatten().tolist()
            del cholesky[2] # this entry is always zero, so don't store it
            values += cholesky
            values.append(bic)
            qmarks = ", ".join(["?"] * len(values))
            try:
                self._cursor.execute(f"INSERT INTO gmm_components VALUES ({qmarks})", values)
            except:
                logger.error(f"Error storing GMM component {i} for kmer {kmer_pos}")
                raise
        self._connection.commit()

    def load_GMM(self, kmer_pos):
        self._cursor.execute("SELECT * FROM gmm_components WHERE kmer_pos = ? ORDER BY component", kmer_pos)
        rows = self._cursor.fetchall()
        if not rows:
            return None # TODO: error?
        model = GaussianMixture(len(rows))
        model.weights_ = np.array([row["weight"] for row in rows])
        model.means_ = np.array([[row["mean_intensity"], row["mean_dwelltime"]] for row in rows])
        model.precisions_cholesky_ = np.array([[[row["precisions_cholesky_00"], row["precisions_cholesky_01"]],
                                                [0.0, row["precisions_cholesky_11"]]] for row in rows])
        # calculate precision matrix from its Cholesky decomposition:
        model.precisions_ = np.array([np.matmul(model.precisions_cholesky_[0],
                                                model.precisions_cholesky_[0].T),
                                      np.matmul(model.precisions_cholesky_[1],
                                                model.precisions_cholesky_[1].T)])
        # calculate covariance matrix from precision matrix:
        model.covariances_ = np.linalg.inv(model.precisions_)
        return model

    def reset_SampComp(self):
        """Reset the database to a state without SampComp results"""
        if not self._connection:
            raise NanocomporeError("Database connection not yet opened")
        # "ALTER table DROP COLUMN" not supported by SQLite version used here
        # (sqlite3.sqlite_version: '3.31.1' - supported from 3.35)
        # can't remove column, so reset values:
        sql = "UPDATE reads SET pass_filter = 0"
        try:
            self._cursor.execute(sql)
        except sqlite3.OperationalError as error:
            # ignore "no such column" error:
            if not error.args[0].startswith("no such column:"):
                raise
        self._cursor.execute("DROP TABLE IF EXISTS gmm_components")
        self._cursor.execute("DROP TABLE IF EXISTS gmm_stats")
        self._cursor.execute("DROP TABLE IF EXISTS kmer_stats")
        self._connection.commit()
