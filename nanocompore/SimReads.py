# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Std lib
from collections import *
import logging
from pkg_resources import resource_filename
import json
import datetime

# Third party
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import pyfaidx
import scipy as sp
import matplotlib.pyplot as pl
from tqdm import tqdm

# Local package
from nanocompore.common import *
from nanocompore import __version__ as package_version
from nanocompore import __name__ as package_name

# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
log_level_dict = {"debug":logging.DEBUG, "info":logging.INFO, "warning":logging.WARNING}

#~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~#
def SimReads (
    fasta_fn,
    outpath="./",
    prefix="reads",
    run_type = "RNA",
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
    mod_extend_context=2,
    min_mod_dist=6,
    pos_rand_seed=42,
    distr_rand_seed=42,
    log_level="info"):
    """
    fasta_fn: Fasta file containing references to use to generate artificial reads
    outpath: Path to the output folder (default "./")
    prefix: text prefix for all the files generated by the function (default "reads")
    run_type: Define the run type model to import (RNA or DNA)
    ref_list: Restrict the references to the listed IDs
    nreads_per_ref: Number of reads to generate per references
    plot: If true, generate an interactive plot of the trace generated
    intensity_mod_loc: value by which to modify the intensity distribution loc value (mode)
    intensity_mod_scale: value by which to modify the intensity distribution scale value (dispersion)
    dwell_mod_loc: value by which to modify the dwell time distribution loc value (mode)
    dwell_mod_scale: value by which to modify the dwell time distribution scale value (mode)
    mod_reads_freq: Frequency of reads to modify
    mod_bases_freq: Frequency of bases to modify in each read (if possible)
    mod_bases_type: Base for which to modify the signal
    mod_extend_context: number of adjacent base affected by the signal modification following an harmonic serries
    min_mod_dist: Minimal distance between to bases to modify
    pos_rand_seed:Define a seed for randon position picking to get a deterministic behaviour
    distr_rand_seed:Define a seed for randon distribution sampling to get a deterministic behaviour
    log_level: Set the log level. Valid values: warning, info, debug
    """

    # Write options in json log file
    kwargs = locals()
    option_d = OrderedDict()
    option_d["package_name"] = package_name
    option_d["package_version"] = package_version
    option_d["function"] = "simulate_reads_from_fasta"
    option_d["timestamp"] = str(datetime.datetime.now())
    for i, j in kwargs.items():
        option_d[i]=j

    # Set logging level
    logger.setLevel(log_level_dict.get(log_level, logging.WARNING))

    if not os.path.exists(outpath):
        logger.debug("Create output dir: {}".format(outpath))
        mkdir (outpath)

    with open ("{}/{}.log".format(outpath, prefix), "w") as log_fp:
        logger.debug("Write log file")
        json.dump (option_d, log_fp, indent=2)

    # Define model depending on run_type
    if run_type == "RNA":
        logger.info("Import RNA model file")
        model_fn = resource_filename("nanocompore", "models/kmers_model_RNA_r9.4_180mv.tsv")
        model_df = pd.read_csv(model_fn, sep="\t", comment="#", index_col=0)
    else:
        raise NanocomporeError("Only RNA is implemented at the moment")

    # Open fasta file and output files
    logger.info("Read Fasta file and simulate corresponding data")

    with pyfaidx.Fasta(fasta_fn) as fasta,\
         open ("{}/{}.tsv".format(outpath, prefix) , "w") as data_fp,\
         open ("{}/{}.tsv.idx".format(outpath, prefix), "w") as idx_fp,\
         open ("{}/{}_pos.tsv".format(outpath, prefix), "w") as pos_fp:

        # Get all reference names if no ref_list
        if not ref_list:
            ref_list = fasta.keys()

        # Write log and index file header and init
        idx_fp.write ("ref_id\tread_id\tbyte_offset\tbyte_len\n")
        pos_fp.write("ref_id\tmodified_positions\n")

        byte_offset = 0
        # Simulate reference per reference
        for ref_num, ref_id in enumerate(tqdm((ref_list), unit=" References", disable=log_level in ("warning", "debug"))):
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
                    pos_rand_seed=pos_rand_seed,
                    distr_rand_seed=distr_rand_seed)

                # Plot traces if required
                if plot:
                    plot_trace (ref_id, intensity_array, dwell_array, mod_pos_list, nreads_mod)

                # Write options used in log file
                pos_fp.write("{}\t{}\n".format(ref_id, array_join(";", mod_pos_list)))

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
    pos_rand_seed=42,
    distr_rand_seed=42):
    """"""

    # Extra parameters if signal modification required
    if mod_bases_freq:
        # Define number of reads to modify and number not to modify
        nreads_mod = int(np.rint(nreads*mod_reads_freq))
        # Define positions to modify base on mod_base_freq and mod_base_type
        mod_pos_list = find_valid_pos_list (ref_seq, mod_bases_type, mod_bases_freq, min_mod_dist, pos_rand_seed)
        # if the modification context has to be extended
        mod_dict = make_mod_dict (intensity_mod_loc, intensity_mod_scale, dwell_mod_loc, dwell_mod_scale, mod_extend_context)
    else:
        mod_pos_list=[]
        nreads_mod=0

    # Create empty arrays to store reads intensity and dwell
    n_kmers = len(ref_seq)-4
    intensity_array = np.empty(shape=(n_kmers, nreads), dtype=np.float)
    dwell_array = np.empty(shape=(n_kmers, nreads), dtype=np.float)

    np.random.seed(distr_rand_seed)
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

def find_valid_pos_list (ref_seq, mod_bases_type, mod_bases_freq, min_mod_dist, pos_rand_seed=42):
    """"""
    pos_list = []
    for i in range(len(ref_seq)-4):
        if ref_seq[i] == mod_bases_type:
            pos_list.append(i)
    n_samples = int(np.rint(len(pos_list)*mod_bases_freq))
    logger.debug ("\tTry to find {} kmers to modify".format(n_samples))

    i = 0
    while True:
        np.random.seed(pos_rand_seed*i)
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