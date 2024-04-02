from nanocompore.common import VAR_ORDER
from nanocompore.common import NanocomporeError

class KmerData():
    def __init__(self, pos, kmer, data, sample_labels, reads, experiment):
        self._kmer = kmer
        self._data = data
        self._pos = pos
        self._sample_labels = sample_labels
        self._reads = reads
        self._experiment = experiment

    ########## Public ##########


    @property
    def pos(self):
        return self._pos


    @property
    def kmer(self):
        return self._kmer


    @property
    def intensity(self):
        return self._data[:, VAR_ORDER.index('trimmean')]


    @property
    def dwell(self):
        return self._data[:, VAR_ORDER.index('dwell')]


    @property
    def sd(self):
        return self._data[:, VAR_ORDER.index('trimsd')]


    @property
    def sample_labels(self):
        return self._sample_labels


    @property
    def reads(self):
        return self._reads


    @property
    def condition_labels(self):
        return list(map(self._experiment.sample_to_condition, self._sample_labels))


    def get_condition_kmer_data(self, condition_label, data_type=''):
        if data_type not in ('intensity', 'dwell'):
            raise NanocomporeError(f"{data_type} is not currently supported")

        return self._data[self._data['condition'] == condition_label][data_type]


    def get_condition_kmer_intensity_data(self, condition_label):
        return self.intensity[[condition_label == cond for cond in self.condition_labels]]


    def get_condition_kmer_dwell_data(self, condition_label):
        return self.dwell[[condition_label == cond for cond in self.condition_labels]]

