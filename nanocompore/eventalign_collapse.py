import sqlite3
import json
import uuid
import sys
import traceback

import multiprocessing as mp

from contextlib import closing
from pathlib import Path

import numpy as np
import statistics

from loguru import logger
from pyfaidx import Fasta

from nanocompore.kmer import KmerData
from nanocompore.database import CREATE_INTERMEDIARY_KMER_DATA_TABLE_QUERY
from nanocompore.database import CREATE_READS_TABLE_QUERY
from nanocompore.database import CREATE_TRANSCRIPTS_TABLE_QUERY
from nanocompore.database import INSERT_INTERMEDIARY_KMER_DATA_QUERY
from nanocompore.database import INSERT_READS_QUERY
from nanocompore.database import INSERT_TRANSCRIPTS_QUERY
from nanocompore.database import CREATE_READS_ID_INDEX_QUERY
from nanocompore.database import CREATE_KMERS_INDEX_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE


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
    def __init__(self, file, fasta_ref, output, kit, nthreads):
        if nthreads < 2:
            raise ValueError("At least 2 threads required.")

        self._file = file
        self._fasta_ref = fasta_ref
        self._ref_lens = {}
        self._kit = kit
        self._nthreads = nthreads
        self._output = output


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

        if self._file:
            self._file_indexer(task_queue)
        else:
            self._stream_indexer(task_queue)

        # add poison pills to kill the workers
        for _ in range(num_workers):
            task_queue.put(None)

        for worker in workers:
            worker.join()

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
                cursor.execute("INSERT INTO kmer_data SELECT * FROM tmp.kmer_data")
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

        tmp_file_name = 'tmp_' + str(uuid.uuid4())
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
                tmp_file_name = 'tmp_' + str(uuid.uuid4())
                tmp_file = open(tmp_file_name, 'w')

                prev_ref_id = ref_id

            tmp_file.write(line)

        if prev_ref_id:
            task_queue.put((prev_ref_id, tmp_file_name))
            tmp_file.close()


    def _read_data(self, file, fstart, fend):
        # The format is:
        # [
        #  ('read_id', [Kmer_pos1, Kmer_pos2, Kmer_pos3...], # read 1
        #  ('read_id', [Kmer_pos1, Kmer_pos2, Kmer_pos3...], # read 2
        # ]
        data = []

        file.seek(fstart)

        logger.info(f"Starting to read lines in range ({fstart}, {fend})")

        prev_read = None
        prev_pos = None
        while True:
            line = file.readline()
            if file.tell() > fend or not line:
                break

            cols = line.split('\t')
            # The position in the eventalign is
            # the 0-based index of the initial
            # base in the k-mer.
            # We want to convert that to 0-based
            # index of the central (most influential)
            # base of the kmer, according to the model.
            # E.g. if we have position 7 in eventalign,
            # then it indicates that the k-mer's first
            # base is the 8th base of the transcript.
            # If RNA002 is used, because the 4th
            # position of the 5mer is the central one,
            # we want to report the 11th base, which
            # in 0-based indexing will have index 10.
            pos = int(cols[1]) + self._kit.center - 1
            # reference kmer
            kmer = cols[2]
            read = cols[3]
            event_measurements = [float(v) for v in cols[13].split(',')]
            event_dwell = float(cols[8])
            model_kmer = cols[9]

            if read != prev_read:
                kmer_data = Kmer(pos,
                                 event_measurements,
                                 event_dwell,
                                 kmer == model_kmer,
                                 kmer)
                # kmer_data = [pos, event_measurements, event_dwell, kmer == model_kmer, kmer]
                data.append((read, [kmer_data]))
                prev_read = read
                prev_pos = pos
            elif pos != prev_pos:
                kmer_data = Kmer(pos,
                                 event_measurements,
                                 event_dwell,
                                 kmer == model_kmer,
                                 kmer)
                # kmer_data = [pos, event_measurements, event_dwell, kmer == model_kmer, kmer]
                data[-1][1].append(kmer_data)
                prev_pos = pos
            else:
                # Get last read, read_data (second element in tuple), last kmer
                kmer_data = data[-1][1][-1]
                kmer_data.measurements.extend(event_measurements)
                kmer_data.dwell += event_dwell
                if kmer_data.valid:
                    kmer_data.valid = kmer == model_kmer

        return data


    def _worker(self, idx, queue, lock, current_transcript_id, current_reads_offset):
        fasta = Fasta(self._fasta_ref)

        tmp_db = self._output + f".{idx}"
        # Make sure there's no old existing database file.
        # If a previous run was abruptly stopped the database
        # may be corrupted, thus breaking the process when
        # we try to open it.
        db_path = Path(self._output)
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
                    kmers = self._process_ref(ref_id, data, ref_len)
                    read_invalid_kmer_ratios = self._get_reads_invalid_kmer_ratio(kmers, ref_len)

                    all_reads = {read
                                 for kmer in kmers
                                 for read in kmer.reads}

                    lock.acquire()
                    transcript_id = current_transcript_id.value
                    current_transcript_id.value += 1
                    reads_offset = current_reads_offset.value
                    current_reads_offset.value += len(all_reads)
                    lock.release()

                    read_ids = {read: i + reads_offset
                                for i, read in enumerate(all_reads)}

                    rows = self._process_rows_for_writing(transcript_id, kmers, read_ids)

                    reads_data = [(read, idx, read_invalid_kmer_ratios[read])
                                  for read, idx in read_ids.items()]

                    cursor.execute('BEGIN')
                    cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                                   (transcript_id, ref_id))
                    cursor.executemany(INSERT_READS_QUERY, reads_data)
                    cursor.executemany(INSERT_INTERMEDIARY_KMER_DATA_QUERY, rows)
                    cursor.execute('COMMIT')
                    logger.info(f"Wrote transcript {transcript_id} {ref_id} to tmp database.")
                except:
                    exception_msg = traceback.format_exc()
                    logger.error(f"Exception for ref {ref_id}: {exception_msg}")
                finally:
                    if not self._file:
                        file.close()
                        Path(filepath).unlink()


    def _process_ref(self, ref_id, data, ref_len):
        nreads = len(data)

        intensity = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        intensity.fill(np.nan)
        sd = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        sd.fill(np.nan)
        dwell = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        dwell.fill(np.nan)
        valid = np.full((nreads, ref_len), False)
        read_ids = []
        kmers = np.repeat('NNNNN', ref_len)

        for row, (read_id, read_data) in enumerate(data):
            read_ids.append(read_id)
            for kmer_data in read_data:
                pos = kmer_data.pos
                # Note: the choice of using statistics' median
                # instead of numpy's is not arbitrary. It seems
                # that statistics.median is faster for small
                # arrays, which is the expected case here.
                median = statistics.median(kmer_data.measurements)
                mad = statistics.median(np.abs(np.array(kmer_data.measurements) - median))
                intensity[row, pos] = median
                sd[row, pos] = mad
                dwell[row, pos] = kmer_data.dwell
                valid[row, pos] = kmer_data.valid
                kmers[pos] = kmer_data.kmer

        kmer_data_list = self._get_kmer_data_list(ref_id,
                                                  ref_len,
                                                  kmers,
                                                  np.array(read_ids),
                                                  intensity,
                                                  sd,
                                                  dwell,
                                                  valid)
        return kmer_data_list


    def _get_kmer_data_list(self,
                            ref_id,
                            ref_len,
                            kmers,
                            read_ids,
                            intensity,
                            sd,
                            dwell,
                            valid):
        kmer_data_list = []
        for pos in range(ref_len):
            valid_reads = ~np.isnan(intensity[:, pos])

            if np.any(valid_reads):
                kmer_data = KmerData(ref_id,
                                     pos,
                                     kmers[pos],
                                     None,
                                     read_ids[valid_reads],
                                     intensity[valid_reads, pos],
                                     sd[valid_reads, pos],
                                     dwell[valid_reads, pos],
                                     valid[valid_reads, pos],
                                     None)
                kmer_data_list.append(kmer_data)

        return kmer_data_list


    def _process_rows_for_writing(self, transcript_id, kmers, read_ids):
        rows = []
        for kmer in kmers:
            kmer_read_ids = np.array([read_ids[read] for read in kmer.reads],
                                     dtype=READ_ID_TYPE)

            rows.append((transcript_id,
                         kmer.pos,
                         kmer.kmer,
                         kmer_read_ids.tobytes(),
                         kmer.intensity.tobytes(),
                         kmer.sd.tobytes(),
                         kmer.dwell.tobytes()))
        return rows


    def _get_reads_invalid_kmer_ratio(self, kmers_data, ref_len):
        read_valid = {}
        for kmer in kmers_data:
            for read, valid in zip(kmer.reads, kmer.valid):
                if read not in read_valid:
                    read_valid[read] = 0
                if valid:
                    read_valid[read] += 1

        return {read: (ref_len - valid)/ref_len
                for read, valid in read_valid.items()}


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
            cursor.execute(CREATE_INTERMEDIARY_KMER_DATA_TABLE_QUERY)
            cursor.execute(CREATE_READS_TABLE_QUERY)
            cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)


    def _create_db_indices(self):
        with closing(self._get_db(self._output)) as conn,\
             closing(conn.cursor()) as cursor:
            cursor.execute(CREATE_READS_ID_INDEX_QUERY)
            cursor.execute(CREATE_KMERS_INDEX_QUERY)

