from sklearn.model_selection import train_test_split

from nanocompore.common import NanocomporeError

class KmerData:
    def __init__(self, pos, kmer, sample_labels, reads, intensity, intensity_std, dwell, valid, experiment):
        self._kmer = kmer
        self._pos = pos
        self._sample_labels = sample_labels
        self._reads = reads
        self._intensity = intensity
        self._intensity_std = intensity_std
        self._dwell = dwell
        self._valid = valid
        self._experiment = experiment


    @property
    def pos(self):
        return self._pos


    @property
    def kmer(self):
        return self._kmer


    @property
    def intensity(self):
        return self._intensity


    @property
    def dwell(self):
        return self._dwell


    @property
    def sd(self):
        return self._intensity_std


    @property
    def sample_labels(self):
        return self._sample_labels


    @property
    def reads(self):
        return self._reads


    @property
    def valid(self):
        return self._valid


    @property
    def condition_labels(self):
        return [self._experiment.sample_to_condition(s)
                for s in self.sample_labels]


    # def get_condition_kmer_data(self, condition_label, data_type=''):
    #     if data_type not in ('intensity', 'dwell'):
    #         raise NanocomporeError(f"{data_type} is not currently supported")

    #     return self._data[self._data['condition'] == condition_label][data_type]


    def get_condition_kmer_intensity_data(self, condition_label):
        return self.intensity[[condition_label == cond for cond in self.condition_labels]]


    def get_condition_kmer_dwell_data(self, condition_label):
        return self.dwell[[condition_label == cond for cond in self.condition_labels]]


    def subsample_reads(self, nreads, random_state=None):
        # Use sklearn to take a stratified sumsample, such
        # that we get an equal number of reads per condition.
        selected, _ = train_test_split(list(range(len(self.intensity))),
                                       train_size=nreads,
                                       stratify=self.condition_labels,
                                       random_state=random_state)
        return KmerData(self.pos,
                        self.kmer,
                        self.sample_labels[selected],
                        self.reads[selected] if self.reads is not None else None,
                        self.intensity[selected],
                        self.sd[selected],
                        self.dwell[selected],
                        self._valid,
                        self._experiment)

