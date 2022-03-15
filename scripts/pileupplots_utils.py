import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import argparse

from collections import defaultdict
from Bio import SeqIO


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vial_fld", help="input folder containing the data for one vial"
    )
    parser.add_argument("--fig_fld", help="output folder where figure are saved")
    parser.add_argument(
        "--show", help="if specified displays the plot", action="store_true"
    )
    return parser


def show_function(do_show):
    """returns a function to either show or close a figure."""
    if do_show:
        return lambda: plt.show()
    else:
        return lambda: plt.close()


def savefig_function(path):
    """returns a function that saves the figure in the right path with
    the given name"""

    def savefunc(name):
        plt.savefig(path / name)

    return savefunc


def load_pileups(data_path):
    """Given the path to the vial_X folder of interest, this function
    loads all the corresponding pileups. Returns them in a dictionary
    whose keys are timepoint labels."""
    # list of "time_*" folders
    timepoint_fld = sorted(list(data_path.glob("time_*")))
    # define regex to extract timepoint label
    m = re.compile("/time_(.+)$")
    # build and return dictionary of pileup matrices
    pileups = {}
    for tpf in timepoint_fld:
        time_id = m.search(str(tpf)).groups()[0]
        tp_file = tpf / "pileup" / "allele_counts.npz"
        pileups[time_id] = np.load(tp_file)["arr_0"]
    return pileups


def load_reference_genome(data_path):
    """Given the path to the vial_X folder of interest, this function
    loads the reference genome from the first timepoint."""
    timepoint_fld = sorted(list(data_path.glob("time_*")))
    fasta_file = timepoint_fld[0] / "ref_genome.fa"
    with open(fasta_file, "r") as f:
        ref_genome = SeqIO.read(f, format="fasta")
    return ref_genome


def coverage(pileup):
    """Returns the sum of nucleotide reads (no gap or N) per position"""
    return pileup.sum(axis=0)[:4, :].sum(axis=0)


def gap_frequency(pileup):
    """returns the frequency of gaps (n. gaps divided by total number of reads, including N)"""
    tot_pileup = pileup.sum(axis=0)
    return tot_pileup[4] / tot_pileup.sum(axis=0)


def ref_genome_to_nucl_id(ref_genome):
    """return an array of indices (0 to 5), one per nucleotide of a
    DNA sequence, according to the order (A,C,G,T,-,N)"""
    nuc_alpha = ["A", "C", "G", "T", "-", "N"]
    nuc_idx = {n: i for i, n in enumerate(nuc_alpha)}
    nseq = np.array(ref_genome)
    return np.array([nuc_idx[x] for x in nseq])


def consensus_frequency(pileup, ref_genome_idxs, count_threshold=5):
    """Returns the consensus frequencies (n. times a nucleotide different from
    the one on the genome is observed). Gaps are discarded"""
    tot_pileup = pileup.sum(axis=0)[:4, :]  # keep only nucleotides, no N or gap
    assert np.all(ref_genome_idxs < 4), "reference sequence includes gaps or Ns"
    n_consensus = np.choose(ref_genome_idxs, tot_pileup)
    consensus_freq = np.zeros_like(n_consensus, dtype=float)
    position_tot = tot_pileup.sum(axis=0)
    mask = position_tot > 0
    consensus_freq[mask] = n_consensus[mask] / position_tot[mask]
    consensus_freq[~mask] = np.nan
    # set reads with less counts than the threshold to NaN
    consensus_freq[position_tot < count_threshold] = np.nan
    return consensus_freq


def color_dict(pileups):
    """Dictionary that assigns a color to every timepoint."""
    cmap = plt.get_cmap("rainbow")
    T = len(pileups)
    colors = {tp: cmap(n / (T - 1)) for n, tp in enumerate(pileups)}
    return colors


def cumulative_histograms(items_dict, ax, colors, plotmeans, **kwargs):
    """plot the cumulative histogram of a dictionary of items"""
    means = {tp: arr.mean() for tp, arr in items_dict.items()}
    for tp, arr in items_dict.items():
        ax.hist(arr, label=tp, histtype="step", color=colors[tp], **kwargs)
        if plotmeans:
            ax.axvline(means[tp], ls=":", color=colors[tp])
