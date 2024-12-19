import pysam
import random
import numpy as np
import pandas as pd
from loguru import logger

from nanocompore.kmer import KmerData
from nanocompore.common import NanocomporeError, Kit
from nanocompore.common import UNCALLED4_MEASUREMENT_TYPE
from nanocompore.common import is_valid_position
from nanocompore.common import get_pos_kmer


class Uncalled4:
    def __init__(self,
                 config,
                 ref_id,
                 seq):
        self._config = config
        self._seq = seq
        self._ref_id = ref_id


    def kmer_data_generator(self):
        """
        A generator that yields one KmerData object
        per position in the transcript.
        """
        kit = self._config.get_kit()
        ref_len = len(self._seq)

        read_ids = []
        sample_labels = []
        condition_labels = []

        read_count = 0
        for sample, _, bam in self._config.get_sample_pod5_bam_data():
            read_count += pysam.AlignmentFile(bam, 'rb').count(reference=self._ref_id)

        # shape is: (reads, positions, vars)
        reads_tensor = np.zeros((read_count, ref_len, 3))

        read_index = 0
        for sample, _, bam in self._config.get_sample_pod5_bam_data():
            for read in pysam.AlignmentFile(bam, 'rb').fetch(self._ref_id):
                if read.is_secondary or read.is_supplementary:
                    continue

                dwell, intensity, sd = self._get_signal(read)

                reads_tensor[read_index, :, 0] = dwell
                reads_tensor[read_index, :, 1] = intensity
                reads_tensor[read_index, :, 2] = sd

                read_ids.append(read.query_name)
                sample_labels.append(sample)
                condition_labels.append(self._config.sample_to_condition(sample))

                read_index += 1

        # Since we may have skipped secondary or supplementary reads,
        # we need to trim the tensor that we preallocated.
        reads_tensor = reads_tensor[:read_index, :, :]

        read_ids = np.array(read_ids)
        sample_labels = np.array(sample_labels)
        condition_labels = np.array(condition_labels)

        # iterate the transcript positions using 1-based indexing
        for pos in range(ref_len):
            # Ignore positions where part of the k-mer is
            # out of the range.
            if not is_valid_position(pos, ref_len, kit):
                continue

            kmer = get_pos_kmer(pos, self._seq, kit)
            pos_data = reads_tensor[:, pos, :]

            # Remove reads that did not have data for the position
            # or represent an unrecoverable skip event (this happens
            # when the skip event is the first one).
            valid_reads = ~np.isnan(pos_data).any(axis=1) & (pos_data[:, 0] != 0)

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
                           None, # We don't have validity data here
                           self._config)


    def _get_signal(self, read):
        """
        Extracts the dwell, intensity and intensity std
        for a given read.
        """
        kit = self._config.get_kit()

        # ur is a list of start, end pairs of aligned segments, e.g:
        # start1, end1, start2, end2
        # ur uses 0-based indexing and [start,end) intervals
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
        ul = self._resolve_skips(ul)

        intensity = []
        sd = []
        dwell = []
        signal_pos = 0
        # We iterate through the segments and get the signal data
        # for each segment in order to merge them in the end.
        for seg_start_indx in range(0, len(read.get_tag('ur')), 2):
            seg_end_indx = seg_start_indx + 1

            # Get the start and end of the region
            # on the reference transcript to which
            # the current signal segment aligns.
            # This is a 0-based [start, end) interval.
            ref_start = ur[seg_start_indx]
            ref_end = ur[seg_end_indx]

            is_first_segment = seg_start_indx == 0
            is_last_segment = seg_end_indx == len(ur) - 1

            # The kmer acts as a sliding window of
            # size > 1. This means that we have fewer
            # measurements than the number of bases in
            # the read. We attribute each measurement to
            # the base that was more influential for the
            # given kmer. For example in RNA002 the most
            # influential base (which we call "center") is
            # the fourth one. This means that the first three
            # bases of the read and the last one would not
            # have measurements (because the kmer cannot
            # "slide out" of the read). Because Uncalled4
            # includes in the reference regions (of the ur tag)
            # all bases that have influenced the measurements
            # of the electrical current we have to trim them.
            if is_first_segment:
                ref_start += kit.center - 1
            if is_last_segment:
                ref_end -= kit.len - kit.center

            current_step = ref_end - ref_start
            next_signal_pos = signal_pos + current_step

            left_pad = ref_start if is_first_segment else 0
            if is_last_segment:
                right_pad = len(self._seq) - ref_end
            else:
                next_start = ur[seg_end_indx + 1]
                right_pad = next_start - ref_end

            intensity.append(self._pad(uc[signal_pos:next_signal_pos], left_pad, right_pad))
            sd.append(self._pad(ud[signal_pos:next_signal_pos], left_pad, right_pad))
            dwell.append(self._pad(ul[signal_pos:next_signal_pos], left_pad, right_pad))
            signal_pos = next_signal_pos

        return np.concatenate(dwell), np.concatenate(intensity), np.concatenate(sd)


    def _pad(self, arr, prefix, suffix):
        return np.concatenate([np.repeat(np.nan, prefix),
                               arr,
                               np.repeat(np.nan, suffix)])


    def _resolve_skips(self, ul):
        resolved = np.empty(len(ul))
        resolved[0] = ul[0]
        for i in range(1, len(ul)):
            resolved[i] = ul[i] if ul[i] != 0 else resolved[i - 1]
        return resolved


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

