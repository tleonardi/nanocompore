"""
Reads an eventalign tsv file produced by
Nanopolish or f5c and creates a collapsed
version of the file that merges together
events for the same position.
The output is written to a compact SQLite
database file.
"""
import multiprocessing as mp
import os
import sqlite3
import statistics
import sys
import threading
import traceback
import uuid

from jaxtyping import Float, Bool
from contextlib import closing
from pathlib import Path

import numpy as np

from loguru import logger
from pyfaidx import Fasta

from nanocompore.database import CREATE_SIGNAL_DATA_TABLE_QUERY
from nanocompore.database import CREATE_READS_TABLE_QUERY
from nanocompore.database import CREATE_TRANSCRIPTS_TABLE_QUERY
from nanocompore.database import INSERT_READS_QUERY
from nanocompore.database import INSERT_SIGNAL_DATA_QUERY
from nanocompore.database import INSERT_TRANSCRIPTS_QUERY
from nanocompore.database import CREATE_READS_ID_INDEX_QUERY
from nanocompore.database import CREATE_SIGNAL_DATA_INDEX_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import MEASUREMENTS_TYPE
from nanocompore.common import monitor_workers


class Kmer:
    def __init__(self, pos, measurements, dwell, valid, kmer):
        self.pos = pos
        self.measurements = measurements
        self.dwell = dwell
        self.valid = valid
        self.kmer = kmer


class EventalignCollapser:
    """
    Parses an eventalign TSV (either given as a
    file or piped to the STDIN of the current
    process), collapses the events for each
    transcript position and stores the relevant
    data to an intermediary SQLite database.
    """
    def __init__(self, file, fasta_ref, output, nthreads, tmp_dir):
        if nthreads < 2:
            raise ValueError("At least 2 threads required.")

        self._file = file
        self._fasta_ref = fasta_ref
        self._ref_lens = {}
        self._nthreads = nthreads
        self._output = output
        self._tmp_dir = tmp_dir


    def __call__(self):
        manager = mp.Manager()
        # nthreads is num_workers + 1, so using this
        # number as a size of the queue means that
        # if the producer is faster than the workers
        # then there'll always be one task ready
        # for processing once a worker gets free.
        # Note: using limited queue is important
        # for the streaming (pipe) case to prevent
        # flooding the file system with big files.
        task_queue = manager.Queue(maxsize=self._nthreads)
        lock = manager.Lock()
        transcript_id = manager.Value('i', 1)
        reads_offset = manager.Value('i', 0)

        num_workers = self._nthreads - 1

        workers = [mp.Process(target=self._worker,
                              args=(i,
                                    task_queue,
                                    lock,
                                    transcript_id,
                                    reads_offset))
                   for i in range(num_workers)]

        for worker in workers:
            worker.start()

        worker_monitor = threading.Thread(target=monitor_workers, args=(workers, ))
        worker_monitor.start()

        if self._file:
            self._file_indexer(task_queue)
        else:
            self._stream_indexer(task_queue)

        # add poison pills to kill the workers
        for _ in range(num_workers):
            task_queue.put(None)

        for worker in workers:
            worker.join()
        worker_monitor.join()

        logger.info("Merging tmp databases.")

        # Make sure there's no old existing database file.
        # If a previous run was abruptly stopped the database
        # may be corrupted, thus breaking the process when
        # we try to open it.
        db_path = Path(self._output)
        if db_path.is_file():
            db_path.unlink()
        self._create_tables(self._output)

        with closing(self._get_db(self._output)) as conn,\
             closing(conn.cursor()) as cursor:
            for i in range(num_workers):
                tmp_db = self._output + f".{i}"
                cursor.execute(f"ATTACH '{tmp_db}' as tmp")
                cursor.execute("BEGIN")
                cursor.execute("INSERT INTO transcripts SELECT * FROM tmp.transcripts")
                cursor.execute("INSERT INTO reads SELECT * FROM tmp.reads")
                cursor.execute("INSERT INTO signal_data SELECT * FROM tmp.signal_data")
                conn.commit()
                cursor.execute("DETACH DATABASE tmp")
                # Delete the tmp db
                Path(tmp_db).unlink()

        logger.info("Creating database indices.")
        self._create_db_indices()


    def _file_indexer(self, task_queue):
        logger.info("The indexer is starting to parse the input eventalign tsv.")
        with open(self._file, 'r') as file:
            # ignore the header
            file.readline()

            prev_ref_id = None
            ref_start = file.tell()
            current_line_pos = None
            while True:
                current_line_pos = file.tell()
                line = file.readline()
                if not line:
                    break

                cols = line.split('\t')
                ref_id = cols[0]

                if not prev_ref_id:
                    prev_ref_id = ref_id

                if ref_id != prev_ref_id:
                    task_queue.put((prev_ref_id, ref_start, current_line_pos))
                    ref_start = current_line_pos
                    prev_ref_id = ref_id
            if prev_ref_id:
                task_queue.put((prev_ref_id, ref_start, current_line_pos))


    def _stream_indexer(self, task_queue):
        # ignore the header
        sys.stdin.readline()

        tmp_file_name = os.path.join(self._tmp_dir,
                                     f'tmp_{str(uuid.uuid4())}')
        tmp_file = open(tmp_file_name, 'w')
        prev_ref_id = None
        while True:
            line = sys.stdin.readline()
            if not line:
                break

            ref_id, _ = line.split('\t', 1)

            if not prev_ref_id:
                prev_ref_id = ref_id

            if ref_id != prev_ref_id:
                task_queue.put((prev_ref_id, tmp_file_name))

                tmp_file.close()
                tmp_file_name = os.path.join(self._tmp_dir,
                                             f'tmp_{str(uuid.uuid4())}')
                tmp_file = open(tmp_file_name, 'w')

                prev_ref_id = ref_id

            tmp_file.write(line)

        if prev_ref_id:
            task_queue.put((prev_ref_id, tmp_file_name))
            tmp_file.close()


    def _read_data(self, file, fstart, fend) -> list[tuple[str, list[Kmer]]]:
        # The format is:
        # [
        #  ('read_id', [Kmer_pos1, Kmer_pos2, Kmer_pos3...], # read 1
        #  ('read_id', [Kmer_pos1, Kmer_pos2, Kmer_pos3...], # read 2
        # ]
        data = []

        file.seek(fstart)

        logger.info(f"Starting to read lines in range [{fstart}, {fend})")

        prev_read = None
        prev_pos = None
        kmer_data = None
        while True:
            line = file.readline()
            if file.tell() > fend or not line:
                break

            cols = line.split('\t')
            # The position in the eventalign is
            # the 0-based index of the initial
            # base in the k-mer.
            pos = int(cols[1])
            # reference kmer
            kmer = cols[2]
            read = cols[3]
            event_measurements = [float(v) for v in cols[13].split(',')]
            event_dwell = float(cols[8])
            model_kmer = cols[9]

            if read != prev_read:
                # If we have a previous kmer, we can finalize it
                # and calculate the median to store only the result
                # and freeing the memory used for all the sampled
                # measurements.
                if kmer_data:
                    # Note: the choice of using statistics' median
                    # instead of numpy's is not arbitrary. It seems
                    # that statistics.median is faster for small
                    # arrays, which is the expected case here.
                    kmer_data.median_intensity = statistics.median(kmer_data.measurements)
                    del kmer_data.measurements
                kmer_data = Kmer(pos,
                                 event_measurements,
                                 event_dwell,
                                 kmer == model_kmer,
                                 kmer)
                data.append((read, [kmer_data]))
                prev_read = read
                prev_pos = pos
            elif pos != prev_pos:
                if kmer_data:
                    kmer_data.median_intensity = statistics.median(kmer_data.measurements)
                    del kmer_data.measurements
                kmer_data = Kmer(pos,
                                 event_measurements,
                                 event_dwell,
                                 kmer == model_kmer,
                                 kmer)
                data[-1][1].append(kmer_data)
                prev_pos = pos
            else:
                # Get last read, read_data (second element in tuple), last kmer
                kmer_data = data[-1][1][-1]
                kmer_data.measurements.extend(event_measurements)
                kmer_data.dwell += event_dwell
                if kmer_data.valid:
                    kmer_data.valid = kmer == model_kmer
        if kmer_data:
            kmer_data.median_intensity = statistics.median(kmer_data.measurements)
            del kmer_data.measurements

        return data


    def _worker(self, idx, queue, lock, current_transcript_id, current_reads_offset):
        fasta = Fasta(self._fasta_ref)

        tmp_db = self._output + f".{idx}"
        # Make sure there's no old existing database file.
        # If a previous run was abruptly stopped the database
        # may be corrupted, thus breaking the process when
        # we try to open it.
        db_path = Path(tmp_db)
        if db_path.is_file():
            db_path.unlink()
        self._create_tables(tmp_db)

        if self._file:
            file = open(self._file, 'r')

        with closing(self._get_db(tmp_db)) as conn,\
             closing(conn.cursor()) as cursor:
            while True:
                task = queue.get()
                if not task:
                    break

                # If we have an eventalign file we can
                # read only the relevant parts with seek.
                # Else, if the eventalign data is streamed
                # through a pipe, we expect the parser to
                # create tmp files for the workers to process.
                if self._file:
                    ref_id, fstart, fend = task
                else:
                    ref_id, filepath = task
                    file = open(filepath, 'r')
                    fstart = 0
                    fend = np.inf

                try:
                    logger.info(f"Starting to read and process data for ref {ref_id}.")
                    data = self._read_data(file, fstart, fend)

                    ref_len = len(str(fasta[ref_id]))
                    intensity, dwell, reads_invalid_ratio = self._process_ref(ref_id, data, ref_len)

                    all_reads = list(reads_invalid_ratio.keys())

                    lock.acquire()
                    transcript_id = current_transcript_id.value
                    current_transcript_id.value += 1
                    reads_offset = current_reads_offset.value
                    current_reads_offset.value += len(all_reads)
                    lock.release()

                    read_ids = {read: i + reads_offset
                                for i, read in enumerate(all_reads)}

                    rows = self._process_rows_for_writing(transcript_id,
                                                          list(read_ids.values()),
                                                          intensity,
                                                          dwell)


                    reads_data = [(read, idx, reads_invalid_ratio[read])
                                  for read, idx in read_ids.items()]

                    cursor.execute('BEGIN')
                    cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                                   (transcript_id, ref_id))
                    cursor.executemany(INSERT_READS_QUERY, reads_data)
                    cursor.executemany(INSERT_SIGNAL_DATA_QUERY, list(rows))
                    cursor.execute('COMMIT')
                    logger.info(f"Wrote transcript {transcript_id} {ref_id} to tmp database.")
                except Exception:
                    exception_msg = traceback.format_exc()
                    logger.error(f"Exception for ref {ref_id}: {exception_msg}")
                finally:
                    if not self._file:
                        file.close()
                        Path(filepath).unlink()


    def _process_ref(
        self,
        ref_id: str,
        data: list[tuple[str, list[Kmer]]],
        ref_len: int
    ) -> tuple[
            Float[np.ndarray, "reads positions"],
            Float[np.ndarray, "reads positions"],
            dict[str, float]]:
        nreads = len(data)

        intensity = np.empty((nreads, ref_len), dtype=MEASUREMENTS_TYPE)
        intensity.fill(np.nan)
        # We don't use sd for now
        # sd = np.empty((nreads, ref_len), dtype=MEASUREMENTS_TYPE)
        # sd.fill(np.nan)
        dwell = np.empty((nreads, ref_len), dtype=MEASUREMENTS_TYPE)
        dwell.fill(np.nan)
        valid = np.full((nreads, ref_len), False)
        read_ids = []
        starts = np.full((nreads,), np.iinfo(np.uint32).max, dtype=np.uint32)
        ends = np.full((nreads,), 0, dtype=np.uint32)

        for row, (read_id, read_data) in enumerate(data):
            read_ids.append(read_id)
            for kmer_data in read_data:
                pos = kmer_data.pos
                intensity[row, pos] = kmer_data.median_intensity
                dwell[row, pos] = kmer_data.dwell
                valid[row, pos] = kmer_data.valid
                starts[row] = min(starts[row], pos)
                ends[row] = max(ends[row], pos)

        # Calculate the ratio of invalid positions for all reads
        alignment_lengths = ends - starts + 1
        invalid_ratios = (alignment_lengths - valid.sum(axis=1))/alignment_lengths
        reads_invalid_ratio = dict(zip(read_ids, invalid_ratios))

        return intensity, dwell, reads_invalid_ratio


    def _process_rows_for_writing(self, transcript_id, read_ids, intensity, dwell):
        for i in range(intensity.shape[0]):
            yield (transcript_id, read_ids[i], intensity[i], dwell[i])


    def _get_db(self, path):
        conn = sqlite3.connect(path, isolation_level=None)
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = OFF')
        # Use largest possible page size of 65kb
        conn.execute('PRAGMA page_size = 65536')
        # Use cache of 2048 x page_size. That's 128MB of cache.
        conn.execute('PRAGMA cache_size = 2000')
        return conn


    def _create_tables(self, path):
        with closing(self._get_db(path)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)
            cursor.execute(CREATE_READS_TABLE_QUERY)
            cursor.execute(CREATE_SIGNAL_DATA_TABLE_QUERY)


    def _create_db_indices(self):
        with closing(self._get_db(self._output)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_READS_ID_INDEX_QUERY)
            cursor.execute(CREATE_SIGNAL_DATA_INDEX_QUERY)

