# -*- coding: utf-8 -*-

import sqlite3
from loguru import logger
from nanocompore.common import NanocomporeError


class DatabaseWrapper(object):

    def __init__(self, db_path):
        self.__db_path = db_path
        self.__connection = None
        self.__cursor = None

    def __enter__(self):
        try:
            logger.debug("Connecting to database")
            self.__connection = sqlite3.connect(self.__db_path)
            self.__connection.row_factory = sqlite3.Row
            self.__cursor = self.__connection.cursor()
        except:
            logger.error("Error connecting to database")
            raise
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.__connection:
            logger.debug("Closing database connection")
            self.__connection.commit()
            self.__connection.close()
            self.__connection = None
            self.__cursor = None

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
            self.cursor.execute("SELECT * FROM samples" + where)
            for row in self.cursor:
                db_samples[row["id"]] = row["name"]
        except Exception:
            logger.error("Error reading sample names from database")
            raise Exception
        for sample in expected_samples: # check that requested samples are in DB
            if sample not in db_samples.values():
                raise NanocomporeError(f"Sample '{sample}' not present in database")
        return db_samples
