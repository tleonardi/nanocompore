import sqlite3 as sql
import numpy as np
import os, collections, sys

class Eventalign_DB():
    def __init__ (self, sample_path):
        self._db_path = sample_path
        self._enterDB()

###############################################################################
#### Private methods ####
###############################################################################
    def _exitDB(self):
        if self._is_connected():
            self._cursor.close()
            self._connection.close()
            sys.stderr.write(f"Connection to {self._db_path} is closed\n")
        else:
            sys.stderr.write(f"{self._db_path} is not connected")

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
            sys.stderr.write(f"connected to {self._db_path}\n")
        else:
            sys.stderr.write(f"{self._db_path} is not a file")
            sys.exit()
    
    def _getTranscript_id(self, transcript):
        tx_querry = f"SELECT t1.id, t1.name FROM transcripts t1 WHERE t1.name = '{transcript}'"
        for row in self._cursor.execute(tx_querry):
            txid = row[0]
            txname = row[1]
        
        return txid
###############################################################################
#### Public methods ####
###############################################################################
    def getData(self, transcript):
        txid = self._getTranscript_id(transcript)

        data_query = f"SELECT k1.readid, k1.position, r1.transcriptid, k1.median, k1.dwell_time FROM kmers k1 LEFT JOIN reads r1 ON r1.id = k1.readid WHERE r1.transcriptid = '{txid}'"

        dwell_time_key = 'dwell_times'
        intensities_key = 'intensities'
        read_ids_key = 'read_ids'
        data = collections.defaultdict(lambda:{dwell_time_key:[], intensities_key:[], read_ids_key:[]})
        for row in self._cursor.execute(data_query):
            read_id = int(row[0])
            pos = int(row[1])
            dwell_time = float(row[4])
            intensity = float(row[3])

            data[pos][read_ids_key].append(read_id)
            data[pos][intensities_key].append(intensity)
            data[pos][dwell_time_key].append(dwell_time)

        for pos in data:
            data[pos][read_ids_key] = np.array(data[pos][read_ids_key])
            data[pos][intensities_key] = np.array(data[pos][intensities_key])
            data[pos][dwell_time_key] = np.array(data[pos][dwell_time_key])
            
        return data

    def closeDB(self):
        self._exitDB()
