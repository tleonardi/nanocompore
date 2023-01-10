#!/usr/bin/env python

import collections, argparse, sys

import numpy as np
import scipy as sp
import pandas as pd

import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.style.use('ggplot')
###############################################################################
def parse_args(args):
    '''
    parses command line arguments
    '''
    parser = argparse.ArgumentParser(description ='helper script that adds a column to the tsv with the dwell time of previous positions')

    parser.add_argument('--collapsed_tsv','-i', type = str, default = '', nargs ='+',
                        help = 'output tsv from nanocompore Eventalign collapse')

    parser.add_argument('--label', '-l', type = str, default = '', nargs ='+',
                        help = 'experimental label for each tsv file')

    parser.add_argument('--output_name','-o', type = str, default ='',
                        help = 'if left blank will use input name as output base name')

    parser.add_argument('--prev_kmer', '-n', type = int, default = 11,
                        help = 'number of kmers backwards to collect dwell from'
                               ' default is 11')

    return parser.parse_args()
###############################################################################

###############################################################################
def parse_tsv(intsv):
    '''
    
    '''
    
    read = ''
    header = ''
    data = []
    with open(intsv, 'r') as infile:
        for line in infile:
            if line.strip().startswith('#'):
                if read:
                    yield read, header, data
                    read = line.strip()
                    header = ''
                    data = []
                else:
                    read = line.strip()
            elif line.strip().startswith('ref_pos'):
                header = line.strip().split('\t')
            else:
                data.append(line.strip().split('\t'))

    yield read, header, data
###############################################################################

###############################################################################
def calc_prev_dwell(intsv, labels, prev_kmer):
    '''
    
    '''
    #refname: pos: label: {'intesnity':[], dwell:[], prev11:[], prev13':[], prev13:[]}
    all_data = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda:
               {'intensity':[], 'dwell':[], 'prev10':[], 'prev11':[], 'prev12':[], 'prev13':[], 'prev13':[], 'prev14':[]})))

    target = 'YNL208W_mRNA'
    target_pos = 352
    off_target_pos = 200
    for file, label in zip(intsv, labels):
        for read, header, data in parse_tsv(file):
            read = read.strip()
            if read != '#':
                ref = read.split()[-1]
                if ref == target:
                    df = pd.DataFrame(data)
                    df.columns = header
                    df['ref_pos'] = df['ref_pos'].astype(int)
                    df['median'] = df['median'].astype(float)
                    df['dwell_time'] = df['dwell_time'].astype(float)
                    max_pos = max(df['ref_pos'])
                    for i in range(target_pos -3, target_pos +4):
                        if len(df.loc[df['ref_pos'] == i, 'ref_pos']) > 0:
                            pos = df.loc[df['ref_pos'] == i, 'ref_pos'].values[0]
                            intensity = df.loc[df['ref_pos'] == i, 'median'].values[0]
                            dwell = df.loc[df['ref_pos'] == i, 'dwell_time'].values[0]
                            if len(df.loc[df['ref_pos'] == i + 10, 'ref_pos']) > 0:
                                #print(target_pos, target_pos+10)
                                #print(df.loc[df['ref_pos'] == target_pos + 10, 'ref_pos'])
                                #print(read)
                                #print(df.loc[df['ref_pos'] == target_pos + 10, 'dwell_time'])
                                prev10 = df.loc[df['ref_pos'] == i + 10, 'dwell_time'].values[0]
                            else:
                                prev10 = 1000.0
                            if len(df.loc[df['ref_pos'] == i + 11, 'ref_pos']) > 0:
                                prev11 = df.loc[df['ref_pos'] == i + 11, 'dwell_time'].values[0]
                            else:
                                prev11 = 1000.0
                            if len(df.loc[df['ref_pos'] == i + 12, 'ref_pos']) > 0:
                                prev12 = df.loc[df['ref_pos'] == i + 12, 'dwell_time'].values[0]
                            else:
                                prev12 = 1000.0
                            if len(df.loc[df['ref_pos'] == i + 13, 'ref_pos']) > 0:
                                prev13 = df.loc[df['ref_pos'] == i + 13, 'dwell_time'].values[0]
                            else:
                                prev13 = 1000.0
                            if len(df.loc[df['ref_pos'] == i + 14, 'ref_pos']) > 0:
                                prev14 = df.loc[df['ref_pos'] == i + 14, 'dwell_time'].values[0]
                            else:
                                prev14 = 1000.0

                            all_data[ref][pos][label]['intensity'].append(intensity)
                            all_data[ref][pos][label]['dwell'].append(dwell)
                            all_data[ref][pos][label]['prev10'].append(prev10)
                            all_data[ref][pos][label]['prev11'].append(prev11)
                            all_data[ref][pos][label]['prev12'].append(prev12)
                            all_data[ref][pos][label]['prev13'].append(prev13)
                            all_data[ref][pos][label]['prev14'].append(prev14)

                    if len(df.loc[df['ref_pos'] == off_target_pos, 'ref_pos']) > 0:
                        pos = pos = df.loc[df['ref_pos'] == off_target_pos, 'ref_pos'].values[0]
                        intensity = df.loc[df['ref_pos'] == off_target_pos, 'median'].values[0]
                        dwell = df.loc[df['ref_pos'] == off_target_pos, 'dwell_time'].values[0]
                        if len(df.loc[df['ref_pos'] == off_target_pos + 10, 'ref_pos']) > 0:
                            prev10 = df.loc[df['ref_pos'] == off_target_pos + 10, 'dwell_time'].values[0]
                        else:
                            prev10 = 1000.0
                        if len(df.loc[df['ref_pos'] == off_target_pos + 11, 'ref_pos']) > 0:
                            prev11 = df.loc[df['ref_pos'] == off_target_pos + 11, 'dwell_time'].values[0]
                        else:
                            prev11 = 1000.0
                        if len(df.loc[df['ref_pos'] == off_target_pos + 12, 'ref_pos']) > 0:
                            prev12 = df.loc[df['ref_pos'] == off_target_pos + 12, 'dwell_time'].values[0]
                        else:
                            prev12 = 1000.0
                        if len(df.loc[df['ref_pos'] == off_target_pos + 13, 'ref_pos']) > 0:
                            prev13 = df.loc[df['ref_pos'] == off_target_pos + 13, 'dwell_time'].values[0]
                        else:
                            prev13 = 1000.0
                        if len(df.loc[df['ref_pos'] == off_target_pos + 14, 'ref_pos']) > 0:
                            prev14 = df.loc[df['ref_pos'] == off_target_pos + 14, 'dwell_time'].values[0]
                        else:
                             prev14 = 1000.0

                        all_data[ref][pos][label]['intensity'].append(intensity)
                        all_data[ref][pos][label]['dwell'].append(dwell)
                        all_data[ref][pos][label]['prev10'].append(prev10)
                        all_data[ref][pos][label]['prev11'].append(prev11)
                        all_data[ref][pos][label]['prev12'].append(prev12)
                        all_data[ref][pos][label]['prev13'].append(prev13)
                        all_data[ref][pos][label]['prev14'].append(prev14)

    return all_data
###############################################################################

###############################################################################
def plot_correlation(data, outfile):
    '''
    
    '''
    color_cycle = ['black', 'red', 'grey', 'gold', 'blue']
    with PdfPages(outfile) as pdf:
        for ref in data:
            for pos in data[ref]:
                #intesnity vs dwell time
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['dwell']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()

                #intesnity vs previous dwell time 10 kmers back
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['prev10']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time 10 pos upstream')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()

                #intesnity vs previous dwell time 11 kmers back
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['prev11']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time 11 pos upstream')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()

                #intesnity vs previous dwell time 12 kmers back
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['prev12']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time 12 pos upstream')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()

                #intesnity vs previous dwell time 13 kmers back
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['prev13']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time 13 pos upstream')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()

                #intesnity vs previous dwell time 14 kmers back
                fig = plt.figure()
                plt.title(f"{ref} {pos}")
                for c, label in enumerate(data[ref][pos]):
                    plt.scatter(data[ref][pos][label]['intensity'],
                                np.log10(data[ref][pos][label]['prev14']),
                                color= color_cycle[c%len(color_cycle)],
                                label = label)

                plt.grid(True, color='grey')
                plt.xlabel('Median ionic current intesnity')
                plt.ylabel('log10 dwell time 14 pos upstream')
                plt.legend()
                pdf.savefig(facecolor='white')
                plt.close()
###############################################################################
def main(args):
    options = parse_args(args)

    all_data = calc_prev_dwell(options.collapsed_tsv, options.label, options.prev_kmer)

    plot_correlation(all_data, options.output_name)
###############################################################################

if (__name__ == "__main__"):
    main(sys.argv)
    raise SystemExit

