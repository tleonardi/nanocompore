"""
Parser for Uncalled4 resquiggling data.

Uncalled4 adds the resquiggling data as additional tags
in the BAM files. This parser reads the BAM file produced
by Uncalled4 and extracts the data from the tags.
"""

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from typing import Dict, Union

import numpy as np
import pysam

from jaxtyping import Float, Int

from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS
from nanocompore.common import Kit


def get_reads(args):
    bam, ref_id, sample = args
    return [(read, sample) for read in bam.fetch(ref_id)]


class Uncalled4:
    def __init__(self,
                 ref_id: str,
                 ref_len: int,
                 bams: Dict[Union[str, int], pysam.AlignmentFile],
                 kit: Kit):
        """
        Initialize an Uncalled4 parser.

        Parameters
        ----------
        ref_id : str
            Reference id for the transcript that will be parsed.
        ref_len : int
            The length of the reference sequence.
        bams : Dict[Union[str, int], pysam.AlignmentFile]
            Dictionary that maps the sample label to the BAM file for the
            sample, opened as a pysam.AlignmentFile.
        kit : Kit
            The chemistry kit used for sequencing the samples.
        """
        self._ref_len = ref_len
        self._ref_id = ref_id
        self._bams = bams
        self._kit = kit


    def get_data(self) -> tuple[Float[np.ndarray, "positions reads vars"],
                                Int[np.ndarray, "reads"],
                                Int[np.ndarray, "reads"]]:
        """
        Returns the signal data for the transcript's reads.

        Returns
        -------
        tuple[Float[np.ndarray, "positions reads vars"],
              Int[np.ndarray, "reads"],
              Int[np.ndarray, "reads"]]
            Tuple with:
            - Tensor with shape (Positions, Reads, Vars)
              containing the signal measurements.
            - 1D array of size <Reads> with read ids.
            - 1D array of size <Reads> with sample labels.
        """
        read_ids = []
        sample_labels = []

        read_count = sum(bam.count(reference=self._ref_id)
                          for bam in self._bams.values())

        # shape is: (positions, reads, vars)
        reads_tensor = np.full((self._ref_len, read_count, 2), np.nan, dtype=np.float32)

        n_samples = len(self._bams)
        read_index = 0
        with ThreadPoolExecutor(max_workers=n_samples) as executor:
            futures = [executor.submit(get_reads, (bam, self._ref_id, sample))
                       for sample, bam in self._bams.items()]
            for future in as_completed(futures):
                sample_reads = future.result()
                for read, sample in sample_reads:
                    if read.is_secondary or read.is_supplementary:
                        continue

                    self._copy_signal_to_tensor(read, reads_tensor, read_index)

                    read_ids.append(read.query_name)
                    sample_labels.append(sample)

                    read_index += 1

        # Since we may have skipped secondary or supplementary reads,
        # we need to trim the tensor that we preallocated.
        reads_tensor = reads_tensor[:, :read_index, :]

        # Uncalled4 will use the minimum value in int16 (i.e. -32768)
        # to encode NaN values for the measurements. We don't want to
        # keep those as valid measurements, because they can lead to
        # false positives.
        reads_tensor = np.where(reads_tensor == -32768, np.nan, reads_tensor)

        if len(read_ids) == 0:
            return reads_tensor, np.array([]), np.array([])

        read_ids = np.array(read_ids)
        sample_labels = np.array(sample_labels)

        return reads_tensor, read_ids, sample_labels


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

        kit = self._kit

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

