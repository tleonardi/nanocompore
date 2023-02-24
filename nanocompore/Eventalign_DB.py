import sqlite3 as sql
import numpy as np
import os, collections, sys

from loguru import logger
from nanocompore.common import *

class Eventalign_DB():
    def __init__ (self, sample_path):
        self._db_path = sample_path
        try:
            self._enterDB()
        except NanocomporeError as e:
            sys.exit()
###############################################################################
#### Private methods ####
###############################################################################
    def _exitDB(self):
        if self._is_connected():
            self._cursor.close()
            self._connection.close()
            logger.trace(f"Closed the connection to {self._db_path}")
        else:
            raise NanocomporeError(f"{self._db_path} is not connected")

    def _is_connected(self):
     try:
        self._connection.cursor()
        return True
     except Exception as ex:
        return False


    def _enterDB(self):
        if os.path.isfile(self._db_path):
            self._connection = sql.connect(self._db_path)
            self._cursor = self._connection.cursor()
            logger.trace(f"Connected to {self._db_path}")
        else:
            raise NanocomporeError(f"{self._db_path} is not a file")
    
    def _getTranscript_id(self, transcript):
        try:
            self._cursor.execute("SELECT id FROM transcripts WHERE name = ?", [transcript])
            row = self._cursor.fetchone()
            if row is not None:
                return row[0]
            else:
                self._cursor.execute("INSERT INTO transcripts (name) VALUES (?)", [transcript])
                self._connection.commit()
                return self._cursor.lastrowid
        except:
            raise NanocomporeError(f"Failed to insert/look up transcript '{transcript}'")

###############################################################################
#### Public methods ####
###############################################################################
    def getData(self, transcript):
        try:
            txid = self._getTranscript_id(transcript)
        except NanocomporeError as e:
            raise

        data_query = f"SELECT k1.readid, k1.position, r1.transcriptid, k1.median, k1.dwell_time, ks.sequence FROM kmers k1 LEFT JOIN reads r1 ON r1.id = k1.readid LEFT JOIN kmer_sequences ks ON k1.sequenceid = ks.id WHERE r1.transcriptid = '{txid}'"

        dwell_time_key = 'dwell_times'
        intensities_key = 'intensities'
        read_ids_key = 'read_ids'
        kmer_seq_key = 'kmer_sequence'
        data = collections.defaultdict(lambda:{dwell_time_key:[], intensities_key:[], read_ids_key:[], kmer_seq_key:[]})
        try:
            for row in self._cursor.execute(data_query):
                read_id = int(row[0])
                pos = int(row[1])
                dwell_time = float(row[4])
                intensity = float(row[3])
                kmer_sequence = row[5]

                data[pos][read_ids_key].append(read_id)
                data[pos][intensities_key].append(intensity)
                data[pos][dwell_time_key].append(dwell_time)
                data[pos][kmer_seq_key].append(kmer_sequence)

            for pos in data:
                data[pos][read_ids_key] = np.array(data[pos][read_ids_key])
                data[pos][intensities_key] = np.array(data[pos][intensities_key])
                data[pos][dwell_time_key] = np.array(data[pos][dwell_time_key])
                data[pos][kmer_seq_key] = collections.Counter(data[pos][kmer_seq_key]).most_common(1)[0][0]
            
            return data
        except:
            self._exitDB()
            raise NanocomporeError(f"Failed to search {self._db_path} for {transcript} data")

    def closeDB(self):
        self._exitDB()
