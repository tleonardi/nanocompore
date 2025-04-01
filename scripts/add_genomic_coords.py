# This should be executed with python >= 3.13 with GIL disabled. E.g.
# $ PYTHON_GIL=0 python3.13t add_genomic_coords.py out_nanocompore_results.tsv gencode41.annotation.gtf

import argparse
import sys
import sqlite3

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import closing
from pathlib import Path

import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='Add genomic positions to Nanocompore results.')
parser.add_argument('input_tsv', help='The results tsv file from Nanocompore.')
parser.add_argument('gtf', help='The annotation GTF.')
parser.add_argument('-t','--threads', help='Number of threads to use (default: 14).', nargs='?', const=14, default=14, type=int)
parser.add_argument('-b','--batch_size', help='Number of rows in each batch (defalut: 3,000,000).', nargs='?', const=3_000_000, default=3_000_000, type=int)
args = parser.parse_args()

INPUT_FN = args.input_tsv
GTF_FN = args.gtf
MAX_ROWS = args.batch_size
THREADS = args.threads

print(f'Will process batches of {MAX_ROWS} using {THREADS} threads.')


class Exon:
    def __init__(self, chrom, strand, start, end):
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end


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

    columns = None
    read = 0

    print('Parsing the GTF annotation.')
    gtf = pd.read_csv(GTF_FN,
                      sep='\t',
                      comment='#',
                      header=None,
                      names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    gtf = gtf[gtf.feature == 'exon']
    gtf['transcript_id'] = gtf.attribute.str.split('; ').apply(lambda x: x[1]).str.split(' ').apply(lambda x: x[1]).str.strip('"')
    all_exons = defaultdict(list)
    for _, row in gtf.iterrows():
        all_exons[row.transcript_id].append(Exon(row.chr, row.strand, row.start, row.end))

    print(f'Got {len(all_exons)} transcripts from the annotation.')


    print('Starting to read CSV and add genomic positions')
    while True:
        if read == 0:
            data = pd.read_csv(INPUT_FN, sep='\t', skiprows=read, nrows=MAX_ROWS)
            columns = data.columns
        else:
            data = pd.read_csv(INPUT_FN, sep='\t', skiprows=read, nrows=MAX_ROWS, names=columns)

        if data.shape[0] == 0:
            break

        print(f'Processing {data.shape[0]} rows.')

        refs = data.ref_id.str.split('|').apply(lambda x: x[0]).unique()

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
    
        if read == 0:
            data.to_csv(result_fn, sep='\t', index=False)
        else:
            data.to_csv(result_fn, sep='\t', index=False, header=False, mode='a')

        read += data.shape[0]

    for w in writers:
        w.join()

