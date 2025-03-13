"""
Parser for Uncalled4 resquiggling data.

Uncalled4 adds the resquiggling data as additional tags
in the BAM files. This parser reads the BAM file produced
by Uncalled4 and extracts the data from the tags.
"""

import pysam
import numpy as np

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from typing import Dict

from nanocompore.kmer import KmerData
from nanocompore.common import UNCALLED4_MEASUREMENT_TYPE
from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS


def get_reads(args):
    bam, ref_id, sample, condition = args
    for read in bam.fetch(ref_id):
        yield read, sample, condition


class Uncalled4:
    def __init__(self,
                 config,
                 ref_id,
                 seq,
                 bams=None):
        """
        Initialize an Uncalled4 parser.

        Parameters
        ----------
        config : nanocompore.config.Config
            A configuration object that contains all the
            parameters for the experiment.
        ref_id : str
            Reference id for the transcript that will be parsed.
        seq : str
            The sequence of the transcript.
        bams : Dict[str, pysam.AlignmentFile]
            Dictionary that maps the sample label to the BAM file for the
            sample, opened as a pysam.AlignmentFile. Default value is None.
            If no value is provided, the files would be opened during the
            initialization of the object. Providing them allows the caller
            to reuse the same opened file for multiple instances.
        """
        self._config = config
        self._seq = seq
        self._ref_id = ref_id
        if bams:
            self._bams = bams
        else:
            self._bams = {sample: pysam.AlignmentFile(sample_def['bam'], 'rb')
                          for _, cond_def in config.get_data().items()
                          for sample, sample_def in cond_def.items()}


    def get_data(self):
        """
        A generator that yields one KmerData object
        per position in the transcript.
        """
        kit = self._config.get_kit()
        ref_len = len(self._seq)

        read_ids = []
        sample_labels = []
        condition_labels = []

        read_counts = defaultdict(lambda: 0)
        for sample, condition, _ in self._config.get_sample_condition_bam_data():
            count = self._bams[sample].count(reference=self._ref_id)
            read_counts[condition] += min(count, self._config.get_downsample_high_coverage())
        read_count = sum(read_counts.values())

        # shape is: (reads, positions, vars)
        reads_tensor = np.zeros((ref_len, read_count, 3))

        n_samples = len(self._config.get_sample_condition_bam_data())
        condition_counts = defaultdict(lambda: 0)
        read_index = 0
        with ThreadPoolExecutor(max_workers=n_samples) as executor:
            futures = [executor.submit(get_reads, (self._bams[sample], self._ref_id, sample, condition))
                       for sample, condition, _ in self._config.get_sample_condition_bam_data()]
            wait(futures, timeout=10)

            zipped_reads = zip(*[future.result() for future in futures])
            for read_group in zipped_reads:
                for read, sample, condition in read_group:
                    if read.is_secondary or read.is_supplementary:
                        continue

                    if condition_counts[condition] >= self._config.get_downsample_high_coverage():
                        continue

                    condition_counts[condition] += 1

                    self._copy_signal_to_tensor(read, reads_tensor, read_index)

                    read_ids.append(read.query_name)
                    sample_labels.append(sample)
                    condition_labels.append(condition)

                    read_index += 1

        # Since we may have skipped secondary or supplementary reads,
        # we need to trim the tensor that we preallocated.
        reads_tensor = reads_tensor[:, :read_index, :]

        read_ids = np.array(read_ids)
        sample_id_mapper = np.vectorize(self._config.get_sample_ids().get)
        sample_ids = sample_id_mapper(sample_labels)
        condition_id_mapper = np.vectorize(self._config.get_condition_ids().get)
        condition_ids = condition_id_mapper(condition_labels)

        return reads_tensor, sample_ids, condition_ids

        # for pos in range(ref_len - kit.len + 1):
        #     kmer = self._seq[pos:pos+kit.len]
        #     pos_data = reads_tensor[:, pos, :]

        #     # Remove reads that did not have data for the position
        #     # or represent an unrecoverable skip event (this happens
        #     # when the skip event is the first one).
        #     valid_reads = ~np.isnan(pos_data[:, 0]) & (pos_data[:, 0] != 0)

        #     if valid_reads.sum() == 0:
        #         continue

        #     pos_data = pos_data[valid_reads, :].astype(UNCALLED4_MEASUREMENT_TYPE)
        #     pos_sample_labels = sample_labels[valid_reads]
        #     pos_read_ids = read_ids[valid_reads]

        #     yield KmerData(self._ref_id,
        #                    pos, # the position is the start of the kmer
        #                    kmer,
        #                    pos_sample_labels,
        #                    pos_read_ids,
        #                    pos_data[:, 1], # intensity
        #                    pos_data[:, 2], # standard dev
        #                    pos_data[:, 0], # dwell
        #                    None, # We don't have validity data here
        #                    self._config)


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
        # We don't use the std for now
        # ud = np.array(read.get_tag('ud')[::-1])
        # ul is a list of dwell times for the kmers
        # The negative length values are paddings used
        # for each alignment segment.
        ul = np.array([d for d in read.get_tag('ul')[::-1] if d >= 0])
        ul = self._resolve_skips(ul)

        segments = [(ur[start_ind], ur[start_ind+1])
                    for start_ind in range(0, len(ur), 2)]

        signal_start = 0
        for ref_start, ref_end in segments:
            # Uncalled4 assigns the measurement to the
            # position of the base that had the highest
            # contribution (so called central position
            # of the kmer, which depending on the pore
            # may not be in the center).
            # The ranges in the ur tag however include
            # all positions that contributed to the
            # measurements. Hence, when we have N
            # measurements, the interval in ur will
            # cover N+K-1 positions (where K is the
            # kmer size for the given pore model).
            # We cut the K-1 positions from the first
            # and the last intervals.
            if ref_start == ur[0]:
                ref_start += kit.len - kit.center
            if ref_end == ur[-1]:
                ref_end -= kit.center - 1
            # Since we want to report the beginning of the
            # kmer instead of the central (most influential)
            # position we shift the signal.
            # Note that the magnitude of this shift is
            # the number of positions between the center and
            # the end of the kmer, because the RNA is read
            # in 3' to 5' direction, so the kmer is inverted
            # with respect to the 5'-3' sequence.
            ref_start -= kit.len - kit.center
            ref_end -= kit.len - kit.center

            segment_size = ref_end - ref_start
            signal_end = signal_start + segment_size

            # Copy signal values to the reads tensor
            tensor[ref_start:ref_end, read_index, DWELL_POS] = ul[signal_start:signal_end]
            tensor[ref_start:ref_end, read_index, INTENSITY_POS] = uc[signal_start:signal_end]
            # We don't use the std for now
            # tensor[read_index, ref_start:ref_end, 2] = ud[signal_start:signal_end]

            signal_start += segment_size


    def _resolve_skips(self, ul):
        resolved = np.empty(len(ul))
        resolved[0] = ul[0]
        for i in range(1, len(ul)):
            resolved[i] = ul[i] if ul[i] != 0 else resolved[i - 1]
        return resolved

