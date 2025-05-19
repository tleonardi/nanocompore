from typing import Optional

import numpy as np
import pod5
import pysam

from jaxtyping import Float
from loguru import logger
from pkg_resources import resource_filename
from remora import io, refine_signal_map, RemoraError

from nanocompore.common import Kit
from nanocompore.common import NanocomporeError


RNA002_LEVELS_FILE = "models/rna002_5mer_levels_v1.txt"
RNA004_LEVELS_FILE = "models/rna004_9mer_levels_v1.txt"
VAR_ORDER = ['dwell', 'trimmean', 'trimsd']

DEFAULT_MAX_READS = 5000


class Remora:
    """
    Resquiggles reads using ONT's Remora.
    """

    def __init__(self,
                 pod5_path: str,
                 bam_path: str,
                 kit: Optional[Kit]=Kit.RNA004,
                 max_reads: Optional[int]=DEFAULT_MAX_READS):
        """

        Parameters
        ----------
        kit : Optional[Kit]
            The sequencing kit used. By default will use RNA004.

        max_reads : Optional[int]
            Maximum number of reads to resquiggle.
            By default will use no more than 5000.

        pod5_path : str
            The path to the pod5 file with the signal data.

        bam_path : str
            The path to the aligned bam file.
        """
        self._kit = kit
        self._max_reads = max_reads
        try:
            self._pod5 = pod5.Reader(pod5_path)
        except Exception:
            raise NanocomporeError(f"Failed to open pod5 file {pod5}")

        try:
            self._bam = pysam.AlignmentFile(bam_path)
        except Exception:
            raise NanocomporeError(f"Failed to open bam file {bam_path}")

        # Remora requires a kmer model file to resquiggle
        # the data (signal to sequence alignment). This is
        # defined as a singal refiner object in the Remora API.
        # Without the signal refiner, it will not resquiggle,
        # but instead merely return the ionic current stream
        try:
            self._sig_map_refiner = self._check_signal_refiner(kit=kit)
        except Exception:
            raise NanocomporeError("Failed to create the signal map refiner. "
                                   "Check that the kmer model table is up-to-date")

        # TODO test Nanocompore accuracy using seconds per kmer
        # or samples per kmer. Remora returns the number of
        # datapoints per kmer, not the number of seconds the
        # kmer persisted in the sensitive region of the nanopore.
        # This function uses the pod5 api to convert the sampling
        # rate of the sequencing (hz) to time per sample (seconds).
        try:
            self._time_per_sample = self._get_time_per_sample()
        except Exception:
            raise NanocomporeError("Failed to check for sampling rate. "
                                   "Likely something wrong with the pod5 file")


    def get_resquiggled_data(
            self,
            ref_id:str,
            ref_seq:str
    ) -> tuple[Float[np.ndarray, "reads positions"],
               Float[np.ndarray, "reads positions"],
               list[str]]:
        """
        Returns the intensity and dwell time data
        for the resquiggled reads.

        Parameters
        ----------
        ref_id : str
            Reference ID of the transcript.
        ref_seq :
            Reference sequence of the transcript.

        Returns
        -------
        tuple[Float[np.ndarray, "reads positions"],
              Float[np.ndarray, "reads positions"],
              list[str]]
            Tuple with:
              - 2D array with shape (reads, positions)
                containing the intensity values.
              - 2D array with shape (reads, positions)
                containing the dwell time values.
              - list with the qname ids of the reads
        """
        ref_region = io.RefRegion(ctg=ref_id,
                                  strand='+',
                                  start=0,
                                  end=len(ref_seq))

        try:
            # Resquiggle the signal and get the summary metrics for
            # each position of the transcript. The result is a dict
            # with metrics:
            # {metric: <numpy appray with shape (reads, positions)>, ...}
            samples_metrics, bam_reads = self._remora_resquiggle(ref_region)
            return (samples_metrics['trimmean'],
                    samples_metrics['dwell'],
                    bam_reads)
        except RemoraError as e:
            if str(e) == "No reads covering region":
                return None
            raise NanocomporeError("failed to resquiggle with Remora") from e
        except Exception as e:
            raise NanocomporeError("failed to resquiggle with Remora") from e


    @property
    def kmer_size(self):
        return self._kmer_size


    def _check_signal_refiner(self, kit):
        level_table = self._kmer_model_selector(kit=kit)
        self._kmer_size = self._kmer_size_detector(level_table=level_table)

        sig_map_refiner = refine_signal_map.SigMapRefiner(
            kmer_model_filename=level_table,
            scale_iters=0,
            algo='dwell_penalty',
            do_rough_rescale=True,
            do_fix_guage=True,
        )
        logger.debug("sig_map_refiner properly opened")
        return sig_map_refiner


    def _get_time_per_sample(self):
        logger.trace("Attempting to calculate time per sample")
        read = next(self._pod5.reads())
        sample_rate = read.run_info.sample_rate
        logger.trace(f"Sampling rate is {sample_rate} Hz")
        time_per_sample = 1.0/sample_rate
        logger.trace(f"Time per sample is {time_per_sample} seconds")
        return time_per_sample


    def _remora_resquiggle(self, ref_region):
        logger.debug(f"Starting to resquiggle data with Remora API for {ref_region.ctg}")
        samples_metrics, all_bam_reads = io.get_ref_reg_samples_metrics(
            ref_region,
            [(self._pod5, self._bam)],
            metric="dwell_trimmean_trimsd",
            sig_map_refiner=self._sig_map_refiner,
            max_reads=self._max_reads,
            # we only support directRNA for now which is sequenced in 3'->5' direction
            reverse_signal=True,
            signal_type='norm',
        )
        logger.debug(f"Data for {ref_region.ctg} resquiggled")
        return samples_metrics[0], [read.qname for read in all_bam_reads[0]]


    def _kmer_model_selector(self, kit):
        if kit == Kit.RNA002:
            level_table = resource_filename('nanocompore', RNA002_LEVELS_FILE)
        elif kit == Kit.RNA004:
            level_table = resource_filename('nanocompore', RNA004_LEVELS_FILE)
        else:
            raise NotImplementedError(f"Kit {kit} not implemented yet.")

        return level_table


    def _kmer_size_detector(self, level_table):
        with open(level_table, 'r') as infile:
            line = infile.readline().strip().split('\t')
            kmer_size = len(line[0])
        logger.debug(f'The kmer size is {kmer_size}')
        return kmer_size

