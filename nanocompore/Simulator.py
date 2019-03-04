# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import logging

# Third party
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import pyfaidx
import scipy as sp
import matplotlib.pyplot as pl

# Local package
from nanocompore.common import *

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
log_level_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#


def simulate_reads_from_fasta (
    fasta_fn,
    model_file,
    output_file,
    ref_list=None,
    nreads_per_ref=100,
    plot=False,
    intensity_mod_loc=0,
    intensity_mod_scale=0,
    dwell_mod_loc=0,
    dwell_mod_scale=0,
    mod_reads_freq=0,
    mod_bases_freq=0,
    mod_bases_type="A",
    mod_extend_context=0,
    min_mod_dist=6,
    rand_seed=42,
    log_level="info"):
    """
    """

    # Saved args with wich the function was called in dict
    kwargs = locals()

    # Set logging level
    logger.setLevel(log_level_dict.get(log_level, logging.WARNING))

    # Import model file
    logger.info("Import model file")
    model_df = pd.read_csv (model_file, sep="\t", index_col=0, comment="#")

    # Open fasta file and output files
    logger.info("Read Fasta file and simulate corresponding data")

    with pyfaidx.Fasta(fasta_fn) as fasta,\
         open (output_file, "w") as data_fp,\
         open (output_file+".idx", "w") as idx_fp,\
         open (output_file+".mdt", "w") as mdt_fp:

        # Get all reference names if no ref_list
        if not ref_list:
            ref_list = fasta.keys()

        # Write metadata header
        for k, v in kwargs.items():
            mdt_fp.write("# {}:{}\n".format(k,v))

        # Write index file header and init
        idx_fp.write ("ref_id\tread_id\tbyte_offset\tbyte_len\n")
        byte_offset = 0

        # Simulate reference per reference
        for ref_num, ref_id in enumerate(ref_list):
            logger.debug("Process reference {}".format(ref_id))
            try:
                ref_seq = str(fasta[ref_id])

                # Simulate data corresponding to the reference
                intensity_array, dwell_array, mod_pos_list, nreads_mod = simulate_ref_mod_context (
                    ref_seq = ref_seq,
                    model_df = model_df,
                    nreads = nreads_per_ref,
                    intensity_mod_loc = intensity_mod_loc,
                    intensity_mod_scale = intensity_mod_scale,
                    dwell_mod_loc = dwell_mod_loc,
                    dwell_mod_scale = dwell_mod_scale,
                    mod_reads_freq = mod_reads_freq,
                    mod_bases_freq = mod_bases_freq,
                    mod_bases_type = mod_bases_type,
                    mod_extend_context = mod_extend_context,
                    min_mod_dist = min_mod_dist,
                    rand_seed = rand_seed)

                # Plot traces if required
                if plot:
                    plot_trace (ref_id, intensity_array, dwell_array, mod_pos_list, nreads_mod)

                # Write modified positions in metadata is any
                if len(mod_pos_list):
                    mdt_fp.write("{}\t{}\n".format(ref_id, array_join(";", mod_pos_list)))

                # Write output in NanopolishComp like files
                for read_num in range(nreads_per_ref):
                    read_str = "#{}_{}\t{}\n".format(ref_num, read_num, ref_id)
                    read_str += "ref_pos\tref_kmer\tdwell_time\tmedian\n"
                    for ref_pos in range(len(ref_seq)-4):
                        read_str += "{}\t{}\t{}\t{}\n".format(
                            ref_pos,
                            ref_seq[ref_pos:ref_pos+5],
                            dwell_array [ref_pos, read_num],
                            intensity_array [ref_pos, read_num])

                    data_fp.write(read_str)
                    idx_fp.write ("{}\t{}_{}\t{}\t{}\n".format(ref_id, ref_num, read_num, byte_offset, len(read_str)-1))
                    byte_offset += len(read_str)

            except KeyError:
                logger.debug ("Reference {} not found in reference fasta file".format(ref_id))

    logger.info("All done")

def plot_trace (ref_id, intensity_array, dwell_array, mod_pos_list, nreads_mod):
    """"""
    with pl.style.context ("ggplot"):
        fig, axes = pl.subplots(2, 1, figsize=(30,10))

        # Plot intensity data
        for i, line in enumerate(intensity_array.T):
            axes[0].plot(line, alpha=(1/len(intensity_array.T))*2, color="red" if i<nreads_mod else "black")
        axes[0].set_title("Median intensity")
        axes[0].set_xlim(0,len(intensity_array))

        # Plot dwell data
        for i, line in enumerate(dwell_array.T):
            axes[1].plot(line, alpha=(1/len(dwell_array.T))*2, color="red" if i<nreads_mod else "black")
        axes[1].set_title("Dwell time")
        axes[1].set_xlim(0,len(dwell_array))

        # Add lines where the signal is modified
        for x in mod_pos_list:
            axes[0].axvline(x, color="grey", linestyle=":")
            axes[1].axvline(x, color="grey", linestyle=":")

        # tweak plot
        fig.suptitle(ref_id, y=1.02, fontsize=18)
        fig.tight_layout()

def simulate_ref_mod_context (
    ref_seq,
    model_df,
    nreads=100,
    intensity_mod_loc=0,
    intensity_mod_scale=0,
    dwell_mod_loc=0,
    dwell_mod_scale=0,
    mod_reads_freq=0,
    mod_bases_freq=0,
    mod_bases_type="A",
    mod_extend_context=0,
    min_mod_dist=6,
    rand_seed=42):
    """"""

    # Extra parameters if signal modification required
    if mod_bases_freq:
        # Define number of reads to modify and number not to modify
        nreads_mod = int(np.rint(nreads*mod_reads_freq))
        # Define positions to modify base on mod_base_freq and mod_base_type
        mod_pos_list = find_valid_pos_list(ref_seq, mod_bases_type, mod_bases_freq, min_mod_dist, rand_seed)
        # if the modification context has to be extended
        mod_dict = make_mod_dict (intensity_mod_loc, intensity_mod_scale, dwell_mod_loc, dwell_mod_scale, mod_extend_context)
    else:
        mod_pos_list=[]
        nreads_mod=0

    # Create empty arrays to store reads intensity and dwell
    n_kmers = len(ref_seq)-4
    intensity_array = np.empty(shape=(n_kmers, nreads), dtype=np.float)
    dwell_array = np.empty(shape=(n_kmers, nreads), dtype=np.float)

    # Fill in arrays with non modified data per position
    for pos in range(n_kmers):
        kmer_seq =  ref_seq[pos:pos+5]
        kmer_model = model_df.loc[kmer_seq]
        intensity_array[pos] = sp.stats.logistic.rvs(loc=kmer_model["model_intensity_loc"], scale=kmer_model["model_intensity_scale"], size=nreads)
        dwell_array[pos] = sp.stats.wald.rvs(loc=kmer_model["model_dwell_loc"], scale=kmer_model["model_dwell_scale"], size=nreads)

    # If modifications are required, edit the values for randomly picked positions + adjacent positions if a context was given
    if mod_bases_freq:
        for pos in mod_pos_list:
            for i in range(-mod_extend_context, mod_extend_context+1):
                pos_extend = pos+i
                if 0 <= pos_extend < n_kmers:
                    kmer_seq = ref_seq[pos_extend:pos_extend+5]
                    kmer_model = model_df.loc[kmer_seq]

                    if intensity_mod_loc or intensity_mod_scale:
                        intensity_array[pos_extend][0:nreads_mod] = sp.stats.logistic.rvs(
                            loc=kmer_model["model_intensity_loc"] + mod_dict["intensity_loc"][i],
                            scale=kmer_model["model_intensity_scale"] + mod_dict["intensity_scale"][i],
                            size=nreads_mod)

                    if dwell_mod_loc or dwell_mod_scale:
                        dwell_array[pos_extend][0:nreads_mod] = sp.stats.wald.rvs(
                            loc=kmer_model["model_dwell_loc"]+ mod_dict["dwell_loc"][i],
                            scale=kmer_model["model_dwell_scale"] + mod_dict["dwell_scale"][i],
                            size=nreads_mod)

    return (intensity_array, dwell_array, mod_pos_list, nreads_mod)

def find_valid_pos_list (ref_seq, mod_bases_type, mod_bases_freq, min_mod_dist, rand_seed=42):
    """"""
    pos_list = []
    for i in range(len(ref_seq)-4):
        if ref_seq[i] == mod_bases_type:
            pos_list.append(i)
    n_samples = int(np.rint(len(pos_list)*mod_bases_freq))
    logger.debug ("\tTry to find {} kmers to modify".format(n_samples))

    i = 0
    while True:
        np.random.seed(rand_seed*i)
        a = np.random.choice(pos_list, n_samples, replace=False)
        a.sort()
        if np.ediff1d(a).min() >= min_mod_dist :
            logger.debug ("\tFound valid combination for {} kmers".format(n_samples))
            logger.debug("\tmodified positions: {}".format(a))
            break
        if i == 1000:
            i = 0
            n_samples -= 1
        i+=1
    return a

def make_mod_dict (intensity_mod_loc, intensity_mod_scale, dwell_mod_loc, dwell_mod_scale, mod_extend_context):
    """Compute a harmonic series per values depending on the context length"""
    pos_list = list (range(-mod_extend_context, mod_extend_context+1))
    d = OrderedDict()
    d["intensity_loc"] = {i: intensity_mod_loc*(1/(abs(i)+1)) for i in pos_list}
    d["intensity_scale"] = {i: intensity_mod_scale*(1/(abs(i)+1)) for i in pos_list}
    d["dwell_loc"] = {i: dwell_mod_loc*(1/(abs(i)+1)) for i in pos_list}
    d["dwell_scale"] = {i: dwell_mod_scale*(1/(abs(i)+1)) for i in pos_list}
    return d

def array_join (sep, array):
    s = ""
    for i in array:
        s+="{}{}".format(i,sep)
    return s[:-len(sep)]
