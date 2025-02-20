# This should be executed with python >= 3.13 with GIL disabled. E.g.
# $ PYTHON_GIL=0 python3.13t add_genomic_coords.py out_nanocompore_results.tsv gencode41.annotation.sqlite

import argparse
import sys
import sqlite3

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import closing
from pathlib import Path
from threading import Thread, Lock

import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='Add genomic positions to Nanocompore results.')
parser.add_argument('input_tsv', help='The results tsv file from Nanocompore.')
parser.add_argument('gtf_db', help='The annotation SQLite DB produced from the GTF.')
parser.add_argument('-t','--threads', help='Number of threads to use (default: 14).', nargs='?', const=14, default=14)
parser.add_argument('-b','--batch_size', help='Number of rows in each batch (defalut: 3,000,000).', nargs='?', const=3_000_000, default=3_000_000)
args = parser.parse_args()

INPUT_FN = args.input_tsv
GTF_FN = args.gtf_db
MAX_ROWS = args.batch_size
THREADS = args.threads

print(MAX_ROWS)


class Exon:
    def __init__(self, chrom, strand, start, end):
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end


def get_exons(refs):
    if not isinstance(refs, set):
        refs = set(refs)
    params = ','.join('?' for _ in refs)
    query = f"""
    SELECT parent, seqid, strand, start, end
    FROM features f
    JOIN (SELECT parent, child
          FROM relations
          WHERE child LIKE "exon_%" AND
                parent IN ({params})) r
    ON f.id = r.child
    """
    exons = defaultdict(list)
    with closing(sqlite3.connect(GTF_FN)) as conn:
        for row in conn.execute(query, tuple(refs)).fetchall():
            # if row[0] in refs:
            exons[row[0]].append(Exon(*row[1:]))
    return exons


def get_transcript_len(exons):
    length = 0
    for exon in exons:
        length += exon.end - exon.start + 1
    return length


def to_genomic_positions(exons):
    if exons[0].strand == '-':
        return np.concatenate([np.arange(e.end, e.start-1, -1) for e in exons])
    return np.concatenate([np.arange(e.start, e.end+1) for e in exons])


def split(m, n):
    k, m = divmod(m, n)
    return [range(i*k+min(i, m), (i+1)*k+min(i+1, m))
            for i in range(n)]


def set_genomic_coords(args):
    r, data, all_exons = args
    current_ref = None
    tx_to_gx = None
    chrom = None
    strand = None
    gx_positions = []
    chroms = []
    strands = []
    for i, row in data.iloc[r].iterrows():
        ref_id = row.ref_id.split('|')[0]
        if ref_id != current_ref:
            current_ref = ref_id
            exons = all_exons[current_ref]
            transcript_len = get_transcript_len(exons)
            # GTF files use 1-based indexing.
            # Hence, we subtract one to convert to 0-based indices.
            tx_to_gx = dict(zip(range(transcript_len),
                                to_genomic_positions(exons)-1))
            chrom = exons[0].chrom
            strand = exons[0].strand
        chroms.append(chrom)
        strands.append(strand)
        gx_positions.append(int(tx_to_gx[row.pos]))
    data.loc[r, 'chr'] = chroms
    data.loc[r, 'strand'] = strands
    data.loc[r, 'genomicPos'] = gx_positions
        

if __name__ == '__main__':
    writers = []
    lock = Lock()

    columns = None
    read = 0
    print('Starting to read CSV and add genomic positions')
    while True:
        if read == 0:
            data = pd.read_csv(INPUT_FN, sep='\t', skiprows=read, nrows=MAX_ROWS)
            columns = data.columns
        else:
            data = pd.read_csv(INPUT_FN, sep='\t', skiprows=read, nrows=MAX_ROWS, names=columns)

        if data.shape[0] == 0:
            break

        print(f'Read {data.shape[0]} rows.')

        print('Getting exons from the annotation...')
        refs = data.ref_id.str.split('|').apply(lambda x: x[0]).unique()
        all_exons = get_exons(refs)
        print(f'Got {len(all_exons)} exons.')

        data['genomicPos'] = pd.Series(dtype='int')
        data['chr'] = pd.Series(dtype='str')
        data['strand'] = pd.Series(dtype='str')
        with ThreadPoolExecutor(max_workers=THREADS) as executor:
            ranges = split(data.shape[0], THREADS)
            for _ in executor.map(set_genomic_coords,
                                  ((r, data, all_exons) for r in ranges)):
                pass
            
        path = Path(INPUT_FN)

        result_fn = "{0}_{2}{1}".format(Path.joinpath(path.parent, path.stem),
                                        path.suffix,
                                        "gx")

        data['genomicPos'] = data['genomicPos'].astype(int)
    
        def write_fn():
            lock.acquire()
            if read == 0:
                data.to_csv(result_fn, sep='\t', index=False)
            else:
                data.to_csv(result_fn, sep='\t', index=False, header=False, mode='a')
            lock.release()

        writer = Thread(target=write_fn)
        writer.start()
        writers.append(writer)

        read += data.shape[0]

    for w in writers:
        w.join()

