import multiprocessing
import os
import queue
import sqlite3
import sys
import time
import traceback

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from concurrent.futures import wait
from contextlib import closing
from multiprocessing.managers import SyncManager

import numpy as np
import numpy.typing as npt
import pysam
import torch

from jaxtyping import Float, Int
from loguru import logger
from pyfaidx import Fasta

from nanocompore.common import Counter
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import SAMPLE_ID_TYPE
from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS
from nanocompore.common import MOTOR_DWELL_POS
from nanocompore.common import MEASUREMENTS_TYPE
from nanocompore.common import UNCALLED4
from nanocompore.common import TranscriptRow
from nanocompore.common import encode_kmer
from nanocompore.common import get_reads_invalid_ratio
from nanocompore.common import get_references_from_bams
from nanocompore.common import get_pos_kmer
from nanocompore.common import log_init_state
from nanocompore.comparisons import TranscriptComparator
from nanocompore.comparisons import retry
from nanocompore.config import Config
from nanocompore.database import ResultsDB
from nanocompore.database import PreprocessingDB
from nanocompore.postprocessing import Postprocessor
from nanocompore.transcript import Transcript
from nanocompore.uncalled4 import Uncalled4


MAX_PROC_ITERATIONS = 500


class RunCmd(object):
    """
    RunCmd reads the resquiggled and preprocessed data
    and performs comparison between the samples of the two
    conditions for each position on every transcript to
    detect the presence of likely modifications.
    """

    def __init__(self, config: Config, random_state: int=42):
        """
        Initialise a `RunCmd` object and generates a white list of
        references with sufficient coverage for subsequent analysis.
        The retuned object can then be called to start the analysis.
        * config
            Config object containing all the parameters for the analysis.
        """
        logger.info("Checking and initialising RunCmd")

        # Save init options in dict for later
        log_init_state(loc=locals())

        self._config = config
        self._workers_with_device = self._get_workers_with_device()
        self._sample_labels = self._config.get_sample_labels()
        self._random_state = random_state


    def __call__(self):
        """
        The main method that will perform the comparative analysis.
        The function will get transcripts for processing from the input
        database, start worker processes to compare them and write the results
        in the output database, and finally will perform the postprocessing
        to export TSVs with the output.
        """
        logger.info("Starting data processing")

        if self._config.get_progress():
            print("Getting the list of transcripts for processing...")
        transcripts = self._get_transcripts_for_processing()
        db_path = ResultsDB(self._config, init_db=False).db_path
        db_exists = os.path.isfile(db_path)
        if db_exists and self._config.get_result_exists_strategy() == 'continue':

            all_transcripts_num = len(transcripts)
            transcripts = self._filter_processed_transcripts(transcripts)
            already_processed = all_transcripts_num - len(transcripts)
            logger.info(f"Found a total of {all_transcripts_num}. " + \
                        f"{already_processed} have already been processed in previous runs. " + \
                        f"Will process only for the remaining {len(transcripts)}.")
        else:
            # Initialize the database
            ResultsDB(self._config, init_db=True)
            logger.info(f"Found a total of {len(transcripts)} transcripts for processing.")

        task_queue = multiprocessing.JoinableQueue()
        for transcript_ref in transcripts:
            task_queue.put((transcript_ref, 0))
        logger.info("All transcripts have been sent to the queue")

        if self._config.get_progress():
            print("Starting to compare the conditions for all transcripts...")
        logger.info(f"Starting {len(self._workers_with_device)} worker processes.")
        db_lock = multiprocessing.Lock()
        sync_lock = multiprocessing.Lock()
        manager = multiprocessing.Manager()
        num_finished = manager.Value('i', 0)

        if self._config.get_resquiggler() == UNCALLED4:
            worker_class = Uncalled4Worker
        else:
            worker_class = GenericWorker
        workers = {worker: worker_class(worker,
                                        task_queue,
                                        db_lock,
                                        sync_lock,
                                        num_finished,
                                        device,
                                        self._config)
                   for worker, device in self._workers_with_device.items()}

        for worker in workers.values():
            worker.start()

        # helper variables for the progress bar
        max_len = 0
        start_time = time.time()

        # Wait for all workers to finish while printing
        # the progress bar if it was requested.
        # We use this slightly weird pattern of waiting
        # for processes instead of a standard loop with
        # unbounded join, because otherwise we'll have
        # to start a separate thread for the progress bar
        # monitoring and that doesn't work well with the
        # multiprocess Manager.
        while len(workers) > 0:
            workers_to_remove = []
            for worker_id, worker in workers.items():
                # Try to join the worker with 1s timeout
                worker.join(1)
                # If the worker has completed we remove it from
                # the worker pool.
                if worker.exitcode is not None:
                    if worker.exitcode == 0:
                        # If we detect a worker that has terminated without an error
                        # it may mean that there were no more jobs in the queue an
                        # the worker closed itself down, or it can be that the worker
                        # did its maximum allowed iterations. If we're in the second
                        # case we want to restart the worker if there are more tasks
                        # to do.
                        with sync_lock:
                            if not task_queue.empty():
                                device = self._workers_with_device[worker_id]
                                logger.info(f"Restarting worker {worker_id} on device {device}.")
                                new_worker = worker_class(worker_id,
                                                          task_queue,
                                                          db_lock,
                                                          sync_lock,
                                                          num_finished,
                                                          device,
                                                          self._config)
                                new_worker.start()
                                workers[worker_id] = new_worker
                            else:
                                workers_to_remove.append(worker_id)
                    else:
                        logger.error("ERROR: A worker encountered an error "
                                     f"(exitcode: {worker.exitcode}). "
                                     "Will terminate all other workers and stop.")
                        for child in multiprocessing.active_children():
                            child.terminate()
                        sys.exit(1)
                if self._config.get_progress():
                    with sync_lock:
                        num_done = num_finished.value
                    self._print_progress_bar(num_done, len(transcripts), max_len, start_time)
            for worker_id in workers_to_remove:
                del workers[worker_id]

        if self._config.get_progress():
            print("\nAll transcripts have been processed.")
            print("Starting postprocessing...")

        postprocessor = Postprocessor(self._config)
        postprocessor()

        if self._config.get_progress():
            print("Done.")


    def _get_transcripts_for_processing(self):
        if self._config.get_resquiggler() == UNCALLED4:
            return get_references_from_bams(self._config)
        else:
            logger.info("Getting references from the input databases.")
            references = set()
            with ThreadPoolExecutor(max_workers=4) as executor:
                futures = [executor.submit(PreprocessingDB(sample_def['db']).get_references_with_data)
                           for condition_def in self._config.get_data().values()
                           for sample, sample_def in condition_def.items()]
                for future in as_completed(futures):
                    references.update(future.result())
            logger.info(f"Found {len(references)} references.")
            return references


    def _filter_processed_transcripts(self, transcripts: set[TranscriptRow]):
        db = ResultsDB(self._config, init_db=False)
        logger.info(f'Preprocessing db: {db.db_path}')
        existing_transcripts = set(db.get_transcripts())
        return [tx for tx in transcripts if tx.ref_id not in existing_transcripts]


    def _get_workers_with_device(self):
        devices = self._config.get_devices()
        result = {}
        if isinstance(devices, dict):
            worker = 0
            for device, procs in devices.items():
                for _ in range(procs):
                    result[worker] = device
                    worker += 1
        else:
            if not isinstance(devices, list):
                devices = [devices]
            for worker in range(self._config.get_nthreads() - 1):
                result[worker] = devices[worker % len(devices)]
        return result


    def _print_progress_bar(self,
                            num_done: int,
                            total_transcripts: int,
                            max_len: int,
                            start_time: float):
        progress_bar_len = 30
        perc_done = num_done / total_transcripts
        done_chars = int(perc_done*progress_bar_len)
        done_bar = '▮' * done_chars
        remaining_bar = '-' * (progress_bar_len - done_chars)
        running_time = time.time() - start_time
        if num_done == 0:
            eta = 'ETA: N/A'
        else:
            total_expected = running_time * total_transcripts / num_done
            remaining_sec = total_expected - running_time
            hours, minutes, _ = self._format_time(remaining_sec)
            eta = f'ETA: {hours:02d}h:{minutes:02d}m remaining'
        hours, minutes, seconds = self._format_time(running_time)
        running = f'⏱  {hours:02d}h:{minutes:02d}m:{seconds:02d}s'
        speed = num_done / (running_time / 60)
        bar = f'{running} |{done_bar}{remaining_bar}| ' + \
              f'{num_done:,}/{total_transcripts:,} ' + \
              f'({perc_done*100:3.1f}%) [{speed:.1f}tx/m] | {eta}'
        if len(bar) > max_len:
            max_len = len(bar)
        if len(bar) < max_len:
            bar = f'{bar}{(max_len - len(bar)) * " "}'
        print(bar, end='\r', flush=True)
        return max_len


    def _format_time(self, seconds: float):
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return int(hours), int(minutes), int(seconds)


class Worker(multiprocessing.Process):
    def __init__(self,
                 worker_id,
                 task_queue,
                 db_lock,
                 sync_lock,
                 num_finished,
                 device,
                 config):
        """
        Abstract parent for different type of worker processes.
        Different workers may read data from different source,
        but then should process it in the same way.

        The worker will continuously take tasks (transcripts)
        from the task queue, process them and write the results
        to the database. Will quit when the task queue is empty.

        Parameters
        ----------
        worker_id : int
            Id of the worker
        task_queue : multiprocessing.JoinableQueue
            Queue with tasks for processing.
        db_lock : multiprocessing.Lock
            Lock to prevent parallel writing to the SQLite DB.
        sync_lock : multiprocessing.Lock
            Lock for synchronisation.
        num_finished : multiprocessing.managers.SyncManager
            Number of processed transcripts.
        device : str
            Which device to use for computation (e.g. cpu or gpu).
        config : nanocompore.Config
            The configuration object.
        """
        super().__init__()
        self._id = worker_id
        self._task_queue = task_queue
        self._db_lock = db_lock
        self._sync_lock = sync_lock
        self._num_finished = num_finished
        self._device = device
        # Note: we're deviating from our convention
        # to use self._config everywhere for the
        # Config object, because _config is used
        # by Process.
        self._conf = config
        self._db_manager = ResultsDB(self._conf, init_db=False)
        self._fasta_fh = Fasta(self._conf.get_fasta_ref())
        self._comparator = TranscriptComparator(config=self._conf)
        self._kit = self._conf.get_kit()
        logger.info(f"Worker {worker_id} started with pid {os.getpid()}.")


    def run(self):
        # When Python grows the heap to find space
        # for memory allocations, after the objects
        # are deleted and freed, the memory may not
        # be returned back to the OS. Hence, the
        # reserved memory stays at the high-water
        # mark level that was previously taken.
        # Since each worker would stumble at an
        # unusally large transcript (in terms of
        # the amount of data) from time to time,
        # the heap of the worker may grow, but it
        # won't shrink after that, leading to
        # a lot of unused head size. When this happens
        # for many workers, the OS or job scheduler
        # may kill the process because it used too
        # much memory.
        # To ameliorate this issue, we can set a
        # fixed number of iterations for the worker
        # and let the main orchestrator process
        # spawn a new worker when it detects the
        # termination.
        for i in range(MAX_PROC_ITERATIONS):
            # It's weird that we have to use a lock here,
            # but for some reason calling get in parallel from
            # multiple processes can cause a queue.Empty
            # exception to be raised even for non-empty queues...
            # https://github.com/python/cpython/issues/87302
            with self._sync_lock:
                try:
                    logger.info(f"Worker {self._id} getting a task from "
                                f"a queue with size {self._task_queue.qsize()}")
                    msg = self._task_queue.get(block=False)
                except queue.Empty:
                    logger.info(f"Worker {self._id} could not find more "
                                "tasks in the queue and will terminate.")
                    break

            transcript_ref, retries = msg
            if retries > 2:
                logger.error(f"Worker {self._id}: Transcript "
                             f"{transcript_ref.ref_id} was retried "
                             "too many times. Won't retry again.")
                continue
            logger.info(f"Worker {self._id} starts processing {transcript_ref.ref_id}")
            transcript = Transcript(id=transcript_ref.id,
                                    ref_id=transcript_ref.ref_id,
                                    ref_seq=str(self._fasta_fh[transcript_ref.ref_id]))

            try:
                data, samples, conditions = self._read_data(transcript)
                # If we don't have valid reads, just continue
                if data.shape[1] == 0:
                    continue
                prepared_data = self._prepare_data(data, samples, conditions)
                data, samples, conditions, positions = prepared_data
                data, positions = self._filter_low_cov_positions(data, positions, conditions)

                transcript, results = self._comparator.compare_transcript(
                        transcript,
                        data,
                        samples,
                        conditions,
                        positions,
                        self._device)

                if results is None:
                    continue

                seq = transcript.seq
                kmers = results.pos.apply(lambda pos: get_pos_kmer(pos, seq, self._kit))
                results['kmer'] = kmers.apply(encode_kmer)
                with self._db_lock:
                    logger.info(f"Worker {self._id}: Saving the results for {transcript.name}")
                    self._db_manager.save_test_results(transcript, results)
            except torch.OutOfMemoryError:
                task = TranscriptRow(transcript.name, transcript.id)
                logger.info(f"Worker {self._id} returns {transcript_ref.ref_id} to " + \
                             "the queue for later processing because it failed with OutOfMemory.")
                with self._sync_lock:
                    self._task_queue.put((task, retries + 1))
            except Exception as e:
                if 'out of memory' in str(e):
                    task = TranscriptRow(transcript.name, transcript.id)
                    logger.info(f"Worker {self._id} returns {transcript_ref.ref_id} to " + \
                                 "the queue for later processing because it failed with OutOfMemory.")
                    with self._sync_lock:
                        self._task_queue.put((task, retries + 1))
                else:
                    msg = traceback.format_exc()
                    logger.error(f"Error in Worker {self._id} for reference {transcript_ref.ref_id}: {msg}")
            finally:
                logger.info(f"Worker {self._id} is finishing a task.")
                with self._sync_lock:
                    self._num_finished.value += 1
                    self._task_queue.task_done()
        else:
            logger.info(f"Worker {self._id} completed a batch of transcripts. Will restart now.")


    def _filter_low_cov_positions(self, data, positions, conditions):
        min_cov = self._conf.get_min_coverage()
        cov_cond_0 = torch.sum(~data[:, conditions == 0, INTENSITY_POS].isnan(), dim=1)
        cov_cond_1 = torch.sum(~data[:, conditions == 1, INTENSITY_POS].isnan(), dim=1)
        valid_positions = torch.logical_and(cov_cond_0 >= min_cov, cov_cond_1 >= min_cov)
        return data[valid_positions], positions[valid_positions]


    def _prepare_data(self,
                      data: Float[np.ndarray, "positions reads vars"],
                      sample_ids: Int[np.ndarray, "reads"],
                      condition_ids: Int[np.ndarray, "reads"],
    ) -> tuple[Float[torch.Tensor, "positions reads vars"],
               Float[torch.Tensor, "reads"],
               Float[torch.Tensor, "reads"],
               Int[torch.Tensor, "positions"]]:
        """
        Prepare data for comparisons.
        This will convert the data to pytorch tensor
        and add the motor dwell time if necessary.

        Parameters
        ----------
        data : Float[np.ndarray, "positions reads vars"]
            Measurements tensor with shape (Positions, Reads, Vars)
        sample_ids : Int[np.ndarray, "reads"]
            1D array of sample ids with size R
        condition_ids : Int[np.ndarray, "reads"]
            1D array of condition ids with size R

        Returns
        -------
        tuple[Float[torch.Tensor, "positions reads vars"],
              Float[torch.Tensor, "reads"],
              Float[torch.Tensor, "reads"],
              Int[torch.Tensor, "positions"]]
            Tuple with:
            - Tensor with shape (Positions, Reads, Vars) containing the
              measurement data.
            - 1D tensor of sample ids with size R
            - 1D tensor of condition ids with size R
            - 1D tensor of reference positions with size P
        """
        motor_offset = self._conf.get_motor_dwell_offset()
        num_positions = data.shape[0]
        if self._conf.get_motor_dwell_offset() > 0:
            data = np.pad(data,
                          ((0, 0), (0, 0), (0, 1)),
                          mode='constant',
                          constant_values=np.nan)
            end = num_positions - motor_offset
            if end > 0:
                data[:end, :, MOTOR_DWELL_POS] = data[motor_offset:, :, DWELL_POS]

        # Mark positions for which we don't have any reads as invalid.
        valid_positions = np.any(~np.isnan(data[:, :, INTENSITY_POS]), axis=1)
        tensor = torch.tensor(data[valid_positions, :, :],
                              dtype=torch.float32,
                              device=self._device)
        samples = torch.tensor(sample_ids,
                               dtype=torch.int16,
                               device=self._device)
        conditions = torch.tensor(condition_ids,
                                  dtype=torch.int16,
                                  device=self._device)
        positions = torch.arange(data.shape[0],
                                 dtype=torch.int16,
                                 device=self._device)[valid_positions]
        return (tensor, samples, conditions, positions)


    def _read_data(
        self,
        transcript: Transcript
    ) -> tuple[Float[np.ndarray, "positions reads vars"],
               Float[np.ndarray, "reads"],
               Float[np.ndarray, "reads"]]:
        """
        Get read data for the given transcript.

        Parameters
        ----------
        transcript : Transcript
            Transcript for which to acqire data.

        Returns
        -------
        tuple[Float[np.ndarray, "positions reads vars"],
              Float[np.ndarray, "reads"],
              Float[np.ndarray, "reads"]]
            Tuple with (measurements tensor, sample ids, condition ids).
            - measurements tensor: shape (Positions, Reads, Vars)
            - sample ids: 1D array with int ids of the samples
            - condition ids: 1D array with int ids of the conditinos (0 or 1).
        """
        raise NotImplementedError("This method should be overriden by specific " + \
                                  "worker process implementations.")


class Uncalled4Worker(Worker):
    def __init__(self,
                 worker_id,
                 task_queue,
                 db_lock,
                 sync_lock,
                 num_finished,
                 device,
                 config):
        """
        Initialize an Uncalled4 worker. It will use
        BAM files produced with Uncalled4's align
        subcommand for its input data.

        The worker will continuously take tasks (transcripts)
        from the task queue, process them and write the results
        to the database. Will quit when the task queue is empty.

        Parameters
        ----------
        worker_id : int
            Id of the worker
        task_queue : multiprocesing.JoinableQueue
            Queue with tasks for processing.
        db_lock : multiprocesing.Lock
            Lock to prevent parallel writing to the SQLite DB.
        sync_lock : multiprocesing.Lock
            Lock for synchronisation.
        num_finished : multiprocesing.managers.SyncManager
            Number of processed transcripts.
        device : str
            Which device to use for computation (e.g. cpu or gpu).
        config : nanocompore.Config
            The configuration object.
        """
        super().__init__(worker_id,
                         task_queue,
                         db_lock,
                         sync_lock,
                         num_finished,
                         device,
                         config)
        data_def = self._conf.get_sample_condition_bam_data()
        self._bams = {sample: pysam.AlignmentFile(bam, 'rb')
                      for sample, _, bam in data_def}


    def _read_data(
        self,
        transcript: Transcript
    ) -> tuple[Float[np.ndarray, "positions reads vars"],
               Float[np.ndarray, "reads"],
               Float[np.ndarray, "reads"]]:
        """
        Get read data for the given transcript.

        Parameters
        ----------
        transcript : Transcript
            Transcript for which to acqire data.

        Returns
        -------
        tuple[Float[np.ndarray, "positions reads vars"],
              Float[np.ndarray, "reads"],
              Float[np.ndarray, "reads"]]
            Tuple with (measurements tensor, sample ids, condition ids).
            - measurements tensor: shape (Positions, Reads, Vars)
            - sample ids: 1D array with int ids of the samples
            - condition ids: 1D array with int ids of the conditinos (0 or 1).
        """
        uncalled4 = Uncalled4(self._conf,
                              transcript.name,
                              transcript.seq,
                              self._bams)
        data, samples, conditions = uncalled4.get_data()

        invalid_ratios = get_reads_invalid_ratio(data[:, :, INTENSITY_POS])
        valid_reads = invalid_ratios < self._conf.get_max_invalid_kmers_freq()
            
        return data[:, valid_reads], samples[valid_reads], conditions[valid_reads]


    def close(self):
        for _, bam in self._bams:
            bam.close()
        super().close()


class GenericWorker(Worker):
    def __init__(self,
                 worker_id,
                 task_queue,
                 db_lock,
                 sync_lock,
                 num_finished,
                 device,
                 config):
        """
        Initialize a generic worker. It will use SQLite databases
        for its input data.

        The worker will continuously take tasks (transcripts)
        from the task queue, process them and write the results
        to the database. Will quit when the task queue is empty.

        Parameters
        ----------
        worker_id : int
            Id of the worker
        task_queue : multiprocesing.JoinableQueue
            Queue with tasks for processing.
        db_lock : multiprocesing.Lock
            Lock to prevent parallel writing to the SQLite DB.
        sync_lock : multiprocesing.Lock
            Lock for synchronisation.
        num_finished : multiprocesing.managers.SyncManager
            Number of processed transcripts.
        device : str
            Which device to use for computation (e.g. cpu or gpu).
        config : nanocompore.Config
            The configuration object.
        """
        super().__init__(worker_id,
                         task_queue,
                         db_lock,
                         sync_lock,
                         num_finished,
                         device,
                         config)
        data_def = self._conf.get_sample_condition_db_data()
        self._dbs = {sample: sqlite3.connect(db)
                     for sample, _, db in data_def}


    def _read_data(
        self,
        transcript: Transcript
    ) -> tuple[Float[np.ndarray, "positions reads vars"],
               Float[np.ndarray, "reads"],
               Float[np.ndarray, "reads"]]:
        """
        Get read data for the given transcript.

        Parameters
        ----------
        transcript : Transcript
            Transcript for which to acqire data.

        Returns
        -------
        tuple[Float[np.ndarray, "positions reads vars"],
              Float[np.ndarray, "reads"],
              Float[np.ndarray, "reads"]]
            Tuple with (measurements tensor, sample ids, condition ids).
            - measurements tensor: shape (Positions, Reads, Vars)
            - sample ids: 1D array with int ids of the samples
            - condition ids: 1D array with int ids of the conditinos (0 or 1).
        """
        intensities = []
        dwells = []
        sample_labels = []
        condition_labels = []
        n_samples = len(self._dbs)
        # Here we will use multi-threaded reading from
        # the input databases (one thread per sample)
        # and then we'll combine the reads into one
        # large tensor.
        with ThreadPoolExecutor(max_workers=n_samples) as executor:
            futures = [executor.submit(get_transcript_signals_from_db,
                                       (self._dbs[sample],
                                        transcript.id,
                                        self._conf.get_max_invalid_kmers_freq(),
                                        self._conf.get_downsample_high_coverage(),
                                        sample,
                                        condition))
                       for sample, condition, _ in self._conf.get_sample_condition_db_data()]
            for future in as_completed(futures):
                signal_data, sample, condition = future.result()
                # intensity will have shape (Reads, Positions)
                intensity = np.array(np.frombuffer(intensity, dtype=MEASUREMENTS_TYPE)
                                     for intensity, _ in signal_data)
                intensities.append(intensity)
                # dwell will have shape (Reads, Positions)
                dwell = np.array(np.frombuffer(dwell, dtype=MEASUREMENTS_TYPE)
                                 for _, dwell in signal_data)
                dwells.append(dwell)
                sample_labels.extend([sample for _ in range(len(signal_data))])
                condition_labels.extend([condition for _ in range(len(signal_data))])

        # The tensor has shape: (Positions, Reads, Vars)
        # Note, intensity and dwell have shape (Reads, Positions)
        # so we transpose them before stacking them.
        tensor = np.stack([np.concatenate(intensities).T,
                           np.concatenate(dwells).T],
                          axis=2)

        sample_id_mapper = np.vectorize(self._conf.get_sample_ids().get)
        sample_ids = sample_id_mapper(sample_labels)
        condition_id_mapper = np.vectorize(self._conf.get_condition_ids().get)
        condition_ids = condition_id_mapper(condition_labels)

        return tensor, sample_ids, condition_ids


    def close(self):
        for _, db in self._dbs:
            db.close()
        super().close()


def get_transcript_signals_from_db(args):
    (connection,
     transcript_id,
     max_invalid_ratio,
     max_reads,
     sample,
     condition) = args
    # Note: we use config.get_downsample_high_coverage to
    # set a limit on the number of reads per condition.
    # However, it will be pointless to take more reads
    # for a given sample than we wull use for the condition.
    # For this reason, we limit the per-sample reads to
    # the same number, to avoid reading and processing
    # large amount of data unnecessarily. When combining
    # the data from different samples of the same conditions
    # some of those reads would be discarded potentially.
    signal_data = PreprocessingDB.get_signal_data(connection,
                                                  transcript_id,
                                                  max_invalid_ratio,
                                                  max_reads)
    return signal_data, sample, condition

