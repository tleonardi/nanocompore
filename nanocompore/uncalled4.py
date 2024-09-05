import pysam
import random
import numpy as np
import pandas as pd
from loguru import logger

from nanocompore.kmer import KmerData
from nanocompore.common import NanocomporeError, Kit
from nanocompore.common import UNCALLED4_MEASUREMENT_TYPE


class Uncalled4:
    def __init__(self,
                 experiment,
                 config,
                 ref_id,
                 seq):
        self._experiment = experiment
        self._config = config
        self._seq = seq
        self._ref_id = ref_id


    def kmer_data_generator(self):
        """
        A generator that yields one KmerData object
        per position in the transcript.
        """
        if self._config.get_kit() == Kit.RNA002:
            kmer_len = 5
        elif self._config.get_kit() == Kit.RNA004:
            kmer_len = 9
        else:
            raise NanocomporeError(f"Kit {self._config.get_kit()} not supported.")
        kmer_radius = (kmer_len - 1) // 2

        ref_len = len(self._seq)

        read_ids = []
        sample_labels = []
        condition_labels = []

        read_count = 0
        for sample, _, bam in self._experiment.get_sample_pod5_bam_data():
            read_count += pysam.AlignmentFile(bam, 'rb').count(reference=self._ref_id)

        # shape is: (reads, positions, vars)
        reads_tensor = np.zeros((read_count, ref_len, 3))

        read_index = 0
        for sample, _, bam in self._experiment.get_sample_pod5_bam_data():
            for read in pysam.AlignmentFile(bam, 'rb').fetch(self._ref_id):
                if read.is_secondary or read.is_supplementary:
                    continue

                dwell, intensity, sd = self._get_signal(read, kmer_radius)

                reads_tensor[read_index, :, 0] = dwell
                reads_tensor[read_index, :, 1] = intensity
                reads_tensor[read_index, :, 2] = sd

                read_ids.append(read.query_name)
                sample_labels.append(sample)
                condition_labels.append(self._experiment.sample_to_condition(sample))

                read_index += 1

        # Since we may have skipped secondary or supplementary reads,
        # we need to trim the tensor that we preallocated.
        reads_tensor = reads_tensor[:read_index, :, :]

        read_ids = np.array(read_ids)
        sample_labels = np.array(sample_labels)
        condition_labels = np.array(condition_labels)

        for pos in range(ref_len):
            if pos < kmer_radius or pos >= ref_len - kmer_radius:
                continue
            # Subtracting one to pos (to both start and end of the interval)
            # to account for the 0-based indexing of the sequence.
            # Adding one to the end of the slice to account for the exclusive nature of the end index.
            # Together, the delta is [-1; 0].
            kmer = self._seq[pos - kmer_radius - 1:pos + kmer_radius]
            pos_data = reads_tensor[:, pos, :]

            # Remove reads that did not have data for the position.
            # non_zero_reads = ~(pos_data == 0).any(axis=1)
            valid_reads = ~np.isnan(pos_data).any(axis=1)

            pos_data = pos_data[valid_reads, :].astype(UNCALLED4_MEASUREMENT_TYPE)
            pos_sample_labels = sample_labels[valid_reads]
            pos_read_ids = read_ids[valid_reads]

            yield KmerData(pos,
                           kmer,
                           pos_sample_labels,
                           pos_read_ids,
                           pos_data[:, 1],
                           pos_data[:, 2],
                           pos_data[:, 0],
                           self._experiment)


    def _get_signal(self, read, kmer_radius):
        """
        Extracts the dwell, intensity and intensity std
        for a given read.
        """
        # ur is a list of start, end pairs of aligned segments, e.g:
        # start1, end1, start2, end2
        ur = read.get_tag('ur')
        # The signal measurements are ordered in
        # 3' to 5' direction so we invert them.
        # uc is the signal's current intensity
        uc = read.get_tag('uc')[::-1]
        # ud is the standard deviation in the current
        ud = read.get_tag('ud')[::-1]
        # ul is a list of dwell times for the kmers
        # The negative length values are paddings used
        # for each alignment segment.
        ul = [d for d in read.get_tag('ul')[::-1] if d >= 0]

        intensity = []
        sd = []
        dwell = []
        current_pos = 0
        # We iterate through the segments and get the signal data
        # for each segment in order to merge them in the end.
        for seg_start_indx in range(0, len(read.get_tag('ur')), 2):
            seg_end_indx = seg_start_indx + 1

            # Inclusive range for the region of the reference
            # sequence that signal aligns to.
            ref_start = ur[seg_start_indx] + kmer_radius
            ref_end = ur[seg_end_indx] - kmer_radius

            next_pos = ref_end - ref_start

            left_pad = ref_start if seg_start_indx == 0 else 0
            if seg_end_indx == len(ur) - 1:
                right_pad = len(self._seq) - ref_end
            else:
                next_start = ur[seg_end_indx + 1] + kmer_radius
                right_pad = next_start - ref_end

            intensity.append(self._pad(uc[current_pos:next_pos], left_pad, right_pad))
            sd.append(self._pad(ud[current_pos:next_pos], left_pad, right_pad))
            dwell.append(self._pad(ul[current_pos:next_pos], left_pad, right_pad))
            next_pos = current_pos

        return np.concatenate(dwell), np.concatenate(intensity), np.concatenate(sd)


    def _pad(self, arr, prefix, suffix):
        return np.concatenate([np.repeat(np.nan, prefix),
                               arr,
                               np.repeat(np.nan, suffix)])


    def _sample(self, labels, limit, seed=42):
        if seed is not None:
            random.seed(seed)

        mask = np.zeros(len(labels), dtype=bool)
        for label in set(labels):
            label_indices = np.where(labels == label)[0]
            random.shuffle(label_indices)
            label_indices = label_indices[:limit]
            for i in label_indices:
                mask[i] = True

        return mask

