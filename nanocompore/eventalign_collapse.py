import sqlite3

from collections import defaultdict

import numpy as np
import statistics

from pyfaidx import Fasta

from nanocompore.kmer import KmerData
from nanocompore.common import DROP_KMER_DATA_TABLE_QUERY
from nanocompore.common import DROP_READS_TABLE_QUERY
from nanocompore.common import DROP_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import CREATE_INTERMEDIARY_KMER_DATA_TABLE_QUERY
from nanocompore.common import CREATE_READS_TABLE_QUERY
from nanocompore.common import CREATE_TRANSCRIPTS_TABLE_QUERY
from nanocompore.common import INSERT_INTERMEDIARY_KMER_DATA_QUERY
from nanocompore.common import INSERT_READS_QUERY
from nanocompore.common import INSERT_TRANSCRIPTS_QUERY
from nanocompore.common import READ_ID_TYPE
from nanocompore.common import EVENTALIGN_MEASUREMENT_TYPE
from nanocompore.common import Indexer


class EventalignCollapser:
    """
    Parses an eventalign TSV (either given as a
    file or piped to the STDIN of the current
    process), collapses the events for each
    transcript position and stores the relevant
    data to an intermediary SQLite database.
    """
    def __init__(self, file, fasta_ref, output, kit):
        if isinstance(file, str):
            self._file = open(file, 'r')
        else:
            self._file = file
        self._fasta_ref = Fasta(fasta_ref)
        self._ref_lens = {}
        self._kit = kit
        self._db_conn = sqlite3.connect(output)
        self._db_cursor = self._db_conn.cursor()
        self._db_cursor.execute(DROP_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(CREATE_INTERMEDIARY_KMER_DATA_TABLE_QUERY)
        self._db_cursor.execute(DROP_READS_TABLE_QUERY)
        self._db_cursor.execute(CREATE_READS_TABLE_QUERY)
        self._db_cursor.execute(DROP_TRANSCRIPTS_TABLE_QUERY)
        self._db_cursor.execute(CREATE_TRANSCRIPTS_TABLE_QUERY)
        self._current_transcript_id = 1
        self._current_read_id = 1

        # transcript -> read -> pos -> data
        self._data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))


    def __call__(self):
        # ignore the header
        self._file.readline()

        prev_ref_id = None
        while True:
            line = self._file.readline()
            if not line:
                break

            cols = line.split('\t')
            ref_id = cols[0]
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

            # Exceptional case for the first line
            if not prev_ref_id:
                prev_ref_id = ref_id

            # If we've come to a line for a new transcript
            # we can process all the accumulated data for
            # the previous one, save it to the database
            # and delete it from memory.
            if ref_id != prev_ref_id:
                self._finish_transcript(prev_ref_id)
                del self._data[prev_ref_id]
                prev_ref_id = ref_id

            # In any case, we store the information
            # from the current line to the appropriate
            # transcript/read/position memory location.
            kmer_data = self._data[ref_id][read][pos]
            if 'measurements' not in kmer_data:
                kmer_data['measurements'] = [event_measurements]
            else:
                kmer_data['measurements'].append(event_measurements)
            if 'dwell' not in kmer_data:
                kmer_data['dwell'] = event_dwell
            else:
                kmer_data['dwell'] += event_dwell
            # If this is the first event we store the current
            # validity status. Else, we update only if it's
            # so far valid. We want to ensure that if one
            # event is invalid, the whole kmer is invalid.
            if 'valid' not in kmer_data or kmer_data['valid']:
                kmer_data['valid'] = kmer == model_kmer
            kmer_data['kmer'] = kmer

        # finish the last transcript
        if prev_ref_id:
            self._finish_transcript(prev_ref_id)


    def _finish_transcript(self, ref_id):
        nreads = len(self._data[ref_id])
        ref_len = len(str(self._fasta_ref[ref_id]))

        intensity = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        intensity.fill(np.nan)
        sd = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        sd.fill(np.nan)
        dwell = np.empty((nreads, ref_len), dtype=EVENTALIGN_MEASUREMENT_TYPE)
        dwell.fill(np.nan)
        valid = np.full((nreads, ref_len), False)
        read_ids = []
        kmers = np.repeat('NNNNN', ref_len)

        for row, read_id in enumerate(self._data[ref_id].keys()):
            read_ids.append(read_id)
            read_data = self._data[ref_id][read_id]
            for pos, pos_data in read_data.items():
                pos_measurements = np.concatenate(pos_data['measurements'])
                # Note: the choice of using statistics' median
                # instead of numpy's is not arbitrary. It seems
                # that statistics.median is faster for small
                # arrays, which is the expected case here.
                median = statistics.median(pos_measurements)
                mad = statistics.median(np.abs(pos_measurements - median))
                intensity[row, pos-1] = median
                sd[row, pos-1] = mad
                dwell[row, pos-1] = pos_data['dwell']
                valid[row, pos-1] = pos_data['valid']
                kmers[pos-1] = pos_data['kmer']

        kmer_data_list = self._get_kmer_data_list(ref_id,
                                                  ref_len,
                                                  kmers,
                                                  np.array(read_ids),
                                                  intensity,
                                                  sd,
                                                  dwell,
                                                  valid)
        self._write_to_db(ref_id, kmer_data_list)


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
                kmer_data = KmerData(pos,
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


    def _write_to_db(self, ref_id, kmer_data_generator):
        self._db_cursor.executemany(
                INSERT_INTERMEDIARY_KMER_DATA_QUERY,
                list(self._process_rows_for_writing(ref_id, kmer_data_generator)))
        self._db_conn.commit()


    def _process_rows_for_writing(self, ref_id, kmers_data):
        ref_len = len(self._fasta_ref[ref_id])

        read_invalid_kmer_ratios = self._get_reads_invalid_kmer_ratio(kmers_data)

        self._db_cursor.execute(INSERT_TRANSCRIPTS_QUERY,
                                (ref_id, self._current_transcript_id))

        read_indexer = Indexer(initial_index=self._current_read_id)
        all_reads = {read
                     for kmer in kmers_data
                     for read in kmer.reads}
        new_mappings = read_indexer.add(all_reads)
        reads_data = [(read, idx, read_invalid_kmer_ratios[read])
                      for read, idx in new_mappings]
        self._db_cursor.executemany(INSERT_READS_QUERY, reads_data)

        for kmer_data in kmers_data:
            kmer_read_ids = read_indexer.get_ids(kmer_data.reads)
            read_ids = np.array(kmer_read_ids, dtype=READ_ID_TYPE)

            yield (self._current_transcript_id,
                   kmer_data.pos,
                   kmer_data.kmer,
                   read_ids.tobytes(),
                   kmer_data.intensity.tobytes(),
                   kmer_data.sd.tobytes(),
                   kmer_data.dwell.tobytes())

        self._current_transcript_id += 1
        self._current_read_id = read_indexer.current_id


    def _get_reads_invalid_kmer_ratio(self, kmers_data):
        read_totals = defaultdict(lambda: 0)
        read_invalid = defaultdict(lambda: 0)
        for kmer in kmers_data:
            for read, valid in zip(kmer.reads, kmer.valid):
                read_totals[read] += 1
                if not valid:
                    read_invalid[read] += 1

        return {read: read_invalid[read]/total
                for read, total in read_totals.items()}

