from collections import Counter

import numpy as np


class KmerData:
    def __init__(self,
                 transcript_id,
                 pos,
                 kmer,
                 sample_labels,
                 reads,
                 intensity,
                 intensity_std,
                 dwell,
                 valid,
                 config):
        self._transcript_id = transcript_id
        self._kmer = kmer
        self._pos = pos
        self._sample_labels = sample_labels
        self._reads = reads
        self._intensity = intensity
        self._intensity_std = intensity_std
        self._dwell = dwell
        self._valid = valid
        self._config = config


    @property
    def transcript_id(self):
        return self._transcript_id


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
    def sample_ids(self):
        sample_to_id = self._config.get_sample_ids()
        return np.vectorize(sample_to_id.get)(self._sample_labels)
    

    @property
    def reads(self):
        return self._reads


    @property
    def valid(self):
        return self._valid


    @property
    def condition_labels(self):
        samp_to_cond = self._config.sample_to_condition()
        return [samp_to_cond[s] for s in self.sample_labels]


    @property
    def condition_ids(self):
        samp_to_cond = self._config.sample_to_condition()
        conditions = np.vectorize(samp_to_cond.get)(self._sample_labels)
        cond_ids = self._config.get_condition_ids()
        return np.vectorize(cond_ids.get)(conditions)


    def get_condition_kmer_intensity_data(self, condition_label):
        return self.intensity[[condition_label == cond for cond in self.condition_labels]]


    def get_condition_kmer_dwell_data(self, condition_label):
        return self.dwell[[condition_label == cond for cond in self.condition_labels]]


    def subsample_reads(self, maximum, random_state=None):
        if all(count <= maximum for count in Counter(self.condition_labels).values()):
            return self

        rand = np.random.default_rng(seed=random_state)
        indices = np.arange(0, len(self.reads))
        all_selected = []
        for condition in set(self.condition_labels):
            condition_mask = np.array(self.condition_labels) == condition
            selected = rand.permutation(indices[condition_mask])[:maximum]
            all_selected.extend(selected)

        mask = np.full(len(self.reads), False)
        for i in all_selected:
            mask[i] = True

        return KmerData(self.transcript_id,
                        self.pos,
                        self.kmer,
                        self.sample_labels[mask],
                        self.reads[mask] if self.reads is not None else None,
                        self.intensity[mask],
                        self.sd[mask],
                        self.dwell[mask],
                        self._valid,
                        self._config)

