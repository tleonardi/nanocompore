"""
Parser for Uncalled4 resquiggling data.

Uncalled4 adds the resquiggling data as additional tags
in the BAM files. This parser reads the BAM file produced
by Uncalled4 and extracts the data from the tags.
"""

import pysam
import numpy as np

from nanocompore.kmer import KmerData
from nanocompore.common import UNCALLED4_MEASUREMENT_TYPE


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
        for sample, _, bam in self._config.get_sample_condition_bam_data():
            read_count += pysam.AlignmentFile(bam, 'rb').count(reference=self._ref_id)

        # shape is: (reads, positions, vars)
        reads_tensor = np.zeros((read_count, ref_len, 3))

        read_index = 0
        for sample, _, bam in self._config.get_sample_condition_bam_data():
            for read in pysam.AlignmentFile(bam, 'rb').fetch(self._ref_id):
                if read.is_secondary or read.is_supplementary:
                    continue

                self._copy_signal_to_tensor(read, reads_tensor, read_index)

                read_ids.append(read.query_name)
                sample_labels.append(sample)
                condition_labels.append(self._config.sample_to_condition()[sample])

                read_index += 1

        # Since we may have skipped secondary or supplementary reads,
        # we need to trim the tensor that we preallocated.
        reads_tensor = reads_tensor[:read_index, :, :]

        read_ids = np.array(read_ids)
        sample_labels = np.array(sample_labels)
        condition_labels = np.array(condition_labels)

        for pos in range(ref_len - kit.len + 1):
            kmer = self._seq[pos:pos+kit.len]
            pos_data = reads_tensor[:, pos, :]

            # Remove reads that did not have data for the position
            # or represent an unrecoverable skip event (this happens
            # when the skip event is the first one).
            valid_reads = ~np.isnan(pos_data[:, 0]) & (pos_data[:, 0] != 0)

            if valid_reads.sum() == 0:
                continue

            pos_data = pos_data[valid_reads, :].astype(UNCALLED4_MEASUREMENT_TYPE)
            pos_sample_labels = sample_labels[valid_reads]
            pos_read_ids = read_ids[valid_reads]

            yield KmerData(self._ref_id,
                           pos, # the position is the start of the kmer
                           kmer,
                           pos_sample_labels,
                           pos_read_ids,
                           pos_data[:, 1], # intensity
                           pos_data[:, 2], # standard dev
                           pos_data[:, 0], # dwell
                           None, # We don't have validity data here
                           self._config)


    def _copy_signal_to_tensor(self, read, tensor, read_index):
        """
        Extracts the dwell, intensity and intensity std
        for a given read.
        """

        # ur gives the reference positions of the measurements.
        # ur is a list of start, end pairs of aligned segments, e.g:
        # start1, end1, start2, end2
        # ur uses 0-based indexing and [start, end) intervals,
        # but it includes the full spans of the evaluated kmers,
        # hence the number of measurements will be less than the
        # positions in the range. Specifically, we would have
        # N-K+1, where N is the reference positions and K is the
        # kmer size.
        ur = read.get_tag('ur')

        kit = self._config.get_kit()

        # The signal measurements are ordered in
        # 3' to 5' direction so we invert them.
        # uc is the signal's current intensity
        uc = np.array(read.get_tag('uc')[::-1])
        # ud is the standard deviation in the current
        ud = np.array(read.get_tag('ud')[::-1])
        # ul is a list of dwell times for the kmers
        # The negative length values are paddings used
        # for each alignment segment.
        ul = np.array([d for d in read.get_tag('ul')[::-1] if d >= 0])
        ul = self._resolve_skips(ul)

        segments = [(ur[start_ind], ur[start_ind+1])
                    for start_ind in range(0, len(ur), 2)]

        signal_start = 0
        for ref_start, ref_end in segments:
            # ref_end is one position after the last
            # one that contributed to a measurement.
            # We subtract kmer-1 to get the position
            # after the starting position of the last
            # included kmer.
            ref_end -= kit.len - 1

            segment_size = ref_end - ref_start
            signal_end = signal_start + segment_size

            # Masked values are assigned a "null" value
            # of -2^16 - 1. We ignore those positions.
            valid = uc != np.iinfo(np.int16).min
            ul = ul[valid]
            uc = uc[valid]
            ud = ud[valid]

            # Copy signal values to the reads tensor
            tensor[read_index, ref_start:ref_end, 0] = ul[signal_start:signal_end]
            tensor[read_index, ref_start:ref_end, 1] = uc[signal_start:signal_end]
            tensor[read_index, ref_start:ref_end, 2] = ud[signal_start:signal_end]

            signal_start += segment_size


    def _resolve_skips(self, ul):
        resolved = np.empty(len(ul))
        resolved[0] = ul[0]
        for i in range(1, len(ul)):
            resolved[i] = ul[i] if ul[i] != 0 else resolved[i - 1]
        return resolved

