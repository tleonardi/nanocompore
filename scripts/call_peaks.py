#!/usr/bin/env python

import argparse, sys
import numpy as np
from scipy.signal import find_peaks


def parse_args(args):
    '''
    parses command line arguments
    '''
    parser = argparse.ArgumentParser(description ='Performs peakcalling on the GMM P-values from the Nanocompore output tsv, filters for significant sites and converts to bed format')

    parser.add_argument('--infile', '-i', type = str, help = 'input Nanocompore output TSV file', required=True)

    parser.add_argument('--outfile', '-o', type = str, help = 'output BED file', required=True)

    parser.add_argument('--kmer_size', '-k', type=int, default=5,
                        help = 'The expected kmer size from nanocompore (default: %(default)s)')

    parser.add_argument('--p_value_threshold', '-p', type=float, default=0.01,
                        help = 'The P-value threshold is used for the height of the peakcalling algorithm (default: %(default)s)')

    parser.add_argument('--LOR_threshold', '-l', type=float, default=0.5,
                        help = 'The Log Odds Ratio (LOR) threshold used to filter for significant sites (default: %(default)s)')

    parser.add_argument('--dynamic_threshold', '-d', default=False, action='store_true',
                        help='Override the P-value threshold and use a transcript specific dynamic threshold for plotting purposes (default: %(default)s)')

    parser.add_argument('--LOR', '-s', default=False, action='store_true',
                        help='Report the Log Odds Ratio (LOR) in the bed instead of the -log10(P-value) (default: %(default)s)')

    return parser.parse_args()


def convert_nanocompore_tsv_to_bed(options):
    infile = options.infile
    outfile = options.outfile
    with open(outfile, 'w') as out:
        for contig, pos, gmm_pvalue, lor in peakcall_nanocompore_transcript_pvalues(
                infile,
                kmer_size=options.kmer_size,
                p_value_threshold=options.p_value_threshold,
                lor_threshold=options.LOR_threshold,
                dynamic=options.dynamic_threshold):
            start = pos
            end = int(pos)+options.kmer_size
            if options.LOR:
                score = lor
                name = "Nanocompore_LOR"
            else:
                score = gmm_pvalue
                name = "Nanocompore_GMM_log10p_value"
            strand = '+'
            out.write(f"{contig}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def peakcall_nanocompore_transcript_pvalues(infile, kmer_size=5, p_value_threshold=0.01, lor_threshold=0.5, dynamic=False):
    for contig, gmm_pvalues, lors, positions in collect_transcript_information(infile):
        # If the p-value is 0, we set it to 1e-323 (the minimum
        # value that doesn't result in -inf) since log10(0) is
        # invalid.
        gmm_pvalues = np.clip(gmm_pvalues, a_min=1e-323, a_max=1)
        gmm_pvalues = -np.log10(gmm_pvalues)
        gmm_pvalues, lors, positions = fill_missing_positions(gmm_pvalues, lors, positions)

        # BED format requires scores between 0 and 1000.
        gmm_pvalue = np.clip(gmm_pvalues, a_min=0, a_max=1000)

        if dynamic:
            height_value = calc_dynamic_threshold(gmm_pvalues)
        else:
            height_value = -np.log10(p_value_threshold)

        peak_idx = find_peaks(gmm_pvalues, height=height_value, distance=kmer_size)

        for idx in peak_idx[0]:
            if gmm_pvalues[idx] >= height_value and abs(lors[idx]) >= lor_threshold:
                yield contig, positions[idx], gmm_pvalues[idx], lors[idx]


def build_values_by_position(gmm_pvalues, lors, positions):
    values_by_position = {}
    for gmm_pvalue, lor, position in zip(gmm_pvalues, lors, positions):
        position = int(position)
        if np.isnan(gmm_pvalue):
            gmm_pvalue = 0.0
        if np.isnan(lor):
            lor = 0.0
        values_by_position[position] = (gmm_pvalue, lor)
    return values_by_position


def fill_missing_positions(gmm_pvalues, lors, positions):
    if not positions:
        return [], [], []

    values_by_position = build_values_by_position(gmm_pvalues, lors, positions)

    min_position = min(values_by_position)
    max_position = max(values_by_position)

    filled_gmm_pvalues = []
    filled_lors = []
    filled_positions = []

    for position in range(min_position, max_position + 1):
        gmm_pvalue, lor = values_by_position.get(position, (0.0, 0.0))

        filled_positions.append(position)
        filled_gmm_pvalues.append(gmm_pvalue)
        filled_lors.append(lor)
    return filled_gmm_pvalues, filled_lors, filled_positions


def calc_dynamic_threshold(pvalues):
    passing_pvalues = []

    for pvalue in pvalues:
        if pvalue >= 2:
            passing_pvalues.append(pvalue)

    if passing_pvalues:
        dynamic_threshold = np.median(sorted(passing_pvalues))
    else:
        dynamic_threshold = 2

    return dynamic_threshold


def collect_transcript_information(infile):
    current_contig = ''
    gmm_pvalues = []
    lors = []
    positions = []
    for contig, gmm_pvalue, lor, pos in read_nanocompore_tsv(infile):
        if current_contig and contig != current_contig:
            yield current_contig, gmm_pvalues, lors, positions
            current_contig = contig
            gmm_pvalues = [gmm_pvalue]
            lors = [lor]
            positions = [pos]
        else:
            current_contig = contig
            gmm_pvalues.append(gmm_pvalue)
            lors.append(lor)
            positions.append(pos)

    yield current_contig, gmm_pvalues, lors, positions


def read_nanocompore_tsv(infile):
    with open(infile, 'r') as tsv:
        line = tsv.readline().strip().split('\t')
        pvalue_indx = line.index('GMM_chi2_pvalue')
        lor_indx = line.index('GMM_LOR')
        for line in tsv:
            line = line.strip().split('\t')
            try:
                gmm_pvalue = float(line[pvalue_indx])
            except:
                gmm_pvalue = float('NaN')

            try:
                lor = float(line[lor_indx])
            except:
                lor = float('NaN')

            contig = line[3].strip()
            pos = line[0].strip()
            yield contig, gmm_pvalue, lor, pos


def main(args):
    #Parse the inputs args/options
    options = parse_args(args)

    convert_nanocompore_tsv_to_bed(options)


if (__name__ == "__main__"):
    main(sys.argv)
    raise SystemExit
