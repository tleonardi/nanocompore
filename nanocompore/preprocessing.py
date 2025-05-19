import multiprocessing as mp
import os
import sqlite3
import threading

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from contextlib import closing
from pathlib import Path
from typing import Any

import pysam
import numpy as np

from jaxtyping import Float, Int
from loguru import logger
from pyfaidx import Fasta

from nanocompore.common import Kit
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import get_reads_invalid_ratio
from nanocompore.common import get_references_from_bam
from nanocompore.common import monitor_workers
from nanocompore.database import INSERT_READS_QUERY
from nanocompore.database import INSERT_TRANSCRIPTS_QUERY
from nanocompore.database import PreprocessingDB
from nanocompore.remora_wrapper import Remora


class Preprocessor:
    """
    The preprocessor takes data coming from
    the chosen resquiggler (described in the
    configuration file used for running the
    experiment) and prepares it in a common
    format in a SQLite DB to be used in the
    subsequest analysis step.
    """
    def __init__(self, fasta_ref: str, output: str, threads: int):
        self._fasta_ref = fasta_ref
        self._output = output
        self._db = PreprocessingDB(output)
        self._db.setup(self._get_metadata())

        self._current_transcript_id = 1
        self._current_read_id = 1

        # The preprocessor will use the main process
        # and will spawn n - 1 worker processes in
        # order to make sure we use the number of
        # processes specified by the user in the config.
        self._worker_processes = threads - 1
        self._references = self._get_references()


    def __call__(self):
        manager = mp.Manager()
        task_queue = manager.Queue()
        lock = manager.Lock()
        current_transcript_id = manager.Value('i', 1)
        current_read_offset = manager.Value('i', 0)
        processed_transcripts = manager.Value('i', 0)

        logger.info(f"Starting {self._worker_processes} worker processes.")

        workers = [mp.Process(target=self._resquiggle,
                              args=(i,
                                    task_queue,
                                    lock,
                                    current_transcript_id,
                                    current_read_offset,
                                    processed_transcripts,
                                    len(self._references)))
                   for i in range(self._worker_processes)]

        for worker in workers:
            worker.start()

        worker_monitor = threading.Thread(target=monitor_workers, args=(workers, ))
        worker_monitor.start()

        for ref_id in self._references:
            task_queue.put(ref_id)

        # add poison pills to kill the workers
        for _ in range(self._worker_processes):
            task_queue.put(None)

        for worker in workers:
            worker.join()
        worker_monitor.join()

        logger.info("Merging tmp databases.")

        # Merge tmp databases
        tmp_dbs = [self._output + f".{i}"
                   for i in range(self._worker_processes)]
        self._db.merge_in_databases(tmp_dbs, merge_reads=True)

        logger.info("Creating database indices.")
        self._db.create_indices()


    def _get_references(self) -> set[str]:
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "resquiggler preprocessor implementations.")


    def _get_resquiggled_data(
            self,
            ref_id: str,
            ref_seq: str
    ) -> tuple[Float[np.ndarray, "reads positions"],
               Float[np.ndarray, "reads positions"],
               list[str]]:
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "resquiggler preprocessor implementations.")


    def _get_metadata(self) -> dict[str, Any]:
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "resquiggler preprocessor implementations.")


    def _resquiggle(self,
                    idx,
                    task_queue,
                    lock,
                    current_transcript_id,
                    current_read_offset,
                    processed_transcripts,
                    num_transcripts):
        fasta_fh = Fasta(self._fasta_ref)

        tmp_db_out = PreprocessingDB(self._output + f".{idx}")
        tmp_db_out.setup(self._get_metadata())

        with closing(tmp_db_out.connect()) as conn:
            while True:
                ref_id = task_queue.get()
                if ref_id is None:
                    break

                ref_seq = str(fasta_fh[ref_id])
                intensity, dwell, reads = self._get_resquiggled_data(ref_id, ref_seq)
                read_invalid_ratios = get_reads_invalid_ratio(intensity.T)

                lock.acquire()
                transcript_id = current_transcript_id.value
                current_transcript_id.value += 1
                read_offset = current_read_offset.value
                current_read_offset.value += len(read_invalid_ratios)
                lock.release()

                reads_with_ratios = zip(reads, read_invalid_ratios)
                offsetted_reads = [(read_id, i + read_offset, invalid_ratio)
                                   for i, (read_id, invalid_ratio) in enumerate(reads_with_ratios)]
                signal_data_rows = self._process_rows_for_writing(
                        transcript_id,
                        [r[1] for r in offsetted_reads], # read internal ids
                        intensity.astype(np.float32),
                        dwell.astype(np.float32))

                with closing(conn.cursor()) as cursor:
                    cursor.execute("begin")
                    cursor.execute(INSERT_TRANSCRIPTS_QUERY, (transcript_id, ref_id))
                    cursor.executemany(INSERT_READS_QUERY, offsetted_reads)
                    cursor.execute("commit")
                PreprocessingDB.write_signal_data_rows(conn, signal_data_rows)

                lock.acquire()
                processed_transcripts.value += 1
                logger.info(f"Data for reference {processed_transcripts.value}/{num_transcripts} "
                            f"{ref_id} has been saved to the tmp database.")
                lock.release()


    def _process_rows_for_writing(self,
                                  transcript_id: int,
                                  read_ids: list[int],
                                  intensity: Float[np.ndarray, "positions"],
                                  dwell: Float[np.ndarray, "positions"]):
        for i in range(intensity.shape[0]):
            yield (transcript_id,
                   read_ids[i],
                   intensity[i].tobytes(),
                   dwell[i].tobytes())


class RemoraPreprocessor(Preprocessor):
    """
    Preprocess data using the Remora resquiggler.
    Remora resquiggling is done during the
    preprocessing (as oposed to Uncalled4 and
    Eventalign where the preprocessor will only
    populate the DB with signal alignment data
    that have already been produced by the resquiggler).
    The signal data would be transfered
    to an SQLite DB to be used by Nanocompore
    for later analysis.
    """
    def __init__(self,
                 fasta_ref: str,
                 pod5: str,
                 bam: str,
                 output: str,
                 kit: Kit,
                 max_reads: int,
                 threads: int):
        """
        Initialize a Remora preprocessor that will
        resquiggle all reads in the provided bam/pod5 pair.

        Parameters
        ----------
        fasta_ref : str
            Path to the fasta reference containing the transcripts.
        pod5 : str
            Path to the pod5 file containing the raw signal data.
        bam : str
            Path to the aligned bam file with the reads to resquiggle.
        output : str
            Path where the output database would be written to.
        kit : Kit
            The sequencing kit used.
        max_reads : int
            Maximum number of reads to resquiggle for each transcript.
        threads : int
            Number of parallel processes to use.
        """
        self._fasta_ref = fasta_ref
        self._pod5 = pod5
        self._bam = bam
        self._kit = kit
        self._max_reads = max_reads
        super().__init__(fasta_ref, output, threads)


    def _get_references(self) -> set[str]:
        return get_references_from_bam(self._bam)


    def _get_resquiggled_data(
            self,
            ref_id,
            ref_seq
    ) -> tuple[Float[np.ndarray, "reads positions"],
               Float[np.ndarray, "reads positions"],
               list[str]]:
        remora = Remora(self._pod5, self._bam, self._kit, self._max_reads)
        return remora.get_resquiggled_data(ref_id, ref_seq)


    def _get_metadata(self) -> dict[str, Any]:
        return {'fasta_ref': self._fasta_ref,
                'pod5': self._pod5,
                'bam': self._bam,
                'kit': self._kit,
                'max_reads': self._max_reads,
                'resquiggler': 'remora'}

