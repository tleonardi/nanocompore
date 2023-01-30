import itertools, Sample
import numpy as np

class Transcript_Data():
    def __init__(self, tx_name, reference_samples, test_samples, tx_length, max_tx_coverage, min_tx_coverage):
        ########## Public fields ##########
        self.name = tx_name
        self.length = tx_length

        ########## Private fields ##########
        self._reference_samples = []
        self._test_samples = []
        self._all_sample_labels = []
        self._min_coverage = min_tx_coverage

        for sample_label in reference_samples:
            self._all_sample_labels.append(sample_label)
            sample_path = reference_samples[sample_label]
            self._reference_samples.append(Sample.Sample('reference', sample_label, sample_path, self.name, max_tx_coverage))

        for sample_label in test_samples:
            self._all_sample_labels.append(sample_label)
            sample_path = test_samples[sample_label]
            self._test_samples.append(Sample.Sample('test', sample_label, sample_path, self.name, max_tx_coverage))

        if len(set(self._all_sample_labels)) == len(self._all_sample_labels):
            self._all_sample_labels = set(self._all_sample_labels)
        else:
            pass
            #[log issue]

        self._valid_positions = self._calcValidPositions(self._min_coverage)

    def enoughTxCoverage(self, min_coverage = ''):
        if min_coverage:
            _min_coverage = min_coverage
        else:
            _min_coverage = self._min_coverage

        if len(self._valid_positions) >= _min_coverage:
            return True
        else:
            return False

    def _calcValidPositions(self, min_coverage):
        intial_positions = set(np.arange(1, self.length+1).tolist())
        invalid_samples = []
        for sample in self._test_samples + self._reference_samples:
            sample_valid_positions = sample.getValidPositions(min_coverage)
            if sample_valid_positions:
                intial_positions = intial_positions.intersection(sample_valid_positions)
            else:
                invalid_samples.append(sample)
        if invalid_samples:
            #log why
            pass
        return intial_positions

    def getValidPositions(self):
        return sorted(self._valid_positions)

    def getReferenceIntensityData(self, pos):
        return self._getDataFromSamples(self._reference_samples, pos, 'intensity')

    def getTestIntensityData(self, pos):
        return self._getDataFromSamples(self._test_samples, pos, 'intensity')

    def _getDataFromSamples(self, sample_list, pos, data_type):
        data = []
        for sample in sample_list:
            if data_type == 'intensity':
                sample_data = sample.getKmerIntensityData(pos)
            elif data_type == 'dwell':
                sample_data = sample.getKmerDwellData(pos)
            elif data_type == 'label':
                sample_data = sample.getKmerLabelData(pos)
            elif data_type == 'read_ids':
                sample_data = sample.getKmerReadIds(pos)
            else:
                sample_data = None
                #log failure and exit

            data.append(sample_data)

        if data_type == 'intensity' or data_type == 'dwell':
            return np.array(data).ravel()
        else:
            return list(itertools.chain(*data))

    def getReferenceDwellData(self, pos):
        return self._getDataFromSamples(self._reference_samples, pos, 'dwell')

    def getTestDwellData(self, pos):
        return self._getDataFromSamples(self._test_samples, pos, 'dwell')

    def getReferenceLabelData(self, pos):
        return self._getDataFromSamples(self._reference_samples, pos, 'label')

    def getTestLabelData(self, pos):
        return self._getDataFromSamples(self._test_samples, pos, 'label')

    def getReferenceReadIdData(self, pos):
        return self._getDataFromSamples(self._reference_samples, pos, 'read_ids')

    def getTestReadIdData(self, pos):
        return self._getDataFromSamples(self._test_samples, pos, 'read_ids')

    def closeAllDbs(self):
        for sample in self._reference_samples + self._test_samples:
            sample.closeDB()