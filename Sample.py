import Data_DB_manager as db_manager
import numpy as np
import random

class Sample():
    def __init__(self, condition_label, sample_label, sample_path, tx_name, max_coverage):
        ############## Public fields ##############
        self.sample_label = sample_label
        self.condition_label = condition_label
        ############## Private fields ##############
        self._sample_path = sample_path
        self._dbManager = db_manager.Data_DB_manager(sample_path)
        self._max_coverage = max_coverage

        self._data = self._dbManager.getData(tx_name)
        self._read_ids_key = 'read_ids'
        self._intensities_key = 'intensities'
        self._dwell_times_key = 'dwell_times'


        if max_coverage != float('inf'):
            self._downSampleHighCoverage()

###############################################################################
#### Private methods ####
###############################################################################
    def _downSampleHighCoverage(self):
        read_ids, down_sample = self._select_reads()

        if down_sample:
            for pos in self._data:
                ids = []
                intensities = []
                dwell_times = []
                for n, id in enumerate(self._data[pos][self._read_ids_key]):
                    if id in read_ids:
                        ids.append(id)
                        intensities.append(self._data[pos][self._intensities_key][n])
                        dwell_times.append(self._data[pos][self._dwell_times_key][n])

                self._data[pos][self._read_ids_key] = np.array(ids)
                self._data[pos][self._intensities_key] = np.array(intensities)
                self._data[pos][self._dwell_times_key] = np.array(dwell_times)               

    def _select_reads(self):
        all_read_ids = set()
        for pos in self._data:
            pos_read_ids = set(self._data[pos][self._read_ids_key])
            all_read_ids = all_read_ids.union(pos_read_ids)

        down_sample = False
        if len(all_read_ids) > self._max_coverage:
            down_sample = True
            reads_to_use = set(random.sample(all_read_ids, self._max_coverage))
        else:
            down_sample = False
            reads_to_use = all_read_ids
            
        return reads_to_use, down_sample
###############################################################################
#### Public methods ####
###############################################################################
    def getValidPositions(self, min_coverage):
        valid_positions = []
        for pos in self._data:
            if len(self._data[pos][self._read_ids_key]) >= min_coverage:
                valid_positions.append(pos)
        return set(valid_positions)

    def getKmerReadIds(self, pos):
        return self._data[pos][self._read_ids_key]

    def getKmerIntensityData(self, pos):
        return self._data[pos][self._intensities_key]

    def getKmerDwellData(self, pos):
        return (self._data[pos][self._dwell_times_key])

    def getKmerLabelData(self, pos):
        label_list = [self.sample_label] * len(self._data[pos][self._read_ids_key])
        return label_list

    def closeDB(self):
        self._dbManager.closeDB()