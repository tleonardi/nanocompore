import nanocompore.Eventalign_DB as Eventalign_DB

class Data_DB_manager():
    def __init__ (self, sample_path):
        self._db = Eventalign_DB.Eventalign_DB(sample_path)
    
    def getData(self, transcript):
        return self._db.getData(transcript)

    def closeDB(self):
        self._db.closeDB()