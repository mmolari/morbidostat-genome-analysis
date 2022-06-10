import numpy as np
import re
import argparse
import pandas as pd
import gzip
import pickle as pkl
import pathlib as pth

from Bio import SeqIO


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vial_fld", help="input folder containing the data for one vial."
    )
    parser.add_argument(
        "--output_fld", help="output folder where the StatsTable is saved."
    )
    parser.add_argument(
        "--verbose", help="print progress statements", action="store_true"
    )
    return parser


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

    def min_timepoint_fld(pathlist):
        ts = []
        for p in pathlist:
            m = re.search(r"/time_(\d+)$", str(p))
            ts.append(int(m.groups()[0]))
        idx = np.argmin(ts)
        return pathlist[idx]

    timepoint_flds = list(data_path.glob("time_*"))
    mintime_fld = min_timepoint_fld(timepoint_flds)
    fasta_file = mintime_fld / "ref_genome.fa"
    with open(fasta_file, "r") as f:
        ref_genome = SeqIO.read(f, format="fasta")
    return ref_genome


def ref_genome_to_nucl_id(ref_genome):
    """return an array of indices (0 to 5), one per nucleotide of a
    DNA sequence, according to the order (A,C,G,T,-,N)"""
    nuc_alpha = ["A", "C", "G", "T", "-", "N"]
    nuc_idx = {n: i for i, n in enumerate(nuc_alpha)}
    nseq = np.array(ref_genome)
    return np.array([nuc_idx[x] for x in nseq])


def consensus_count(pileup, ref_genome):
    """Returns the number of consensus nucleotides (n. times a nucleotide is
    the one on the reference genome) and the total number of nucleotides. Gaps are discarded"""
    ref_genome_idxs = ref_genome_to_nucl_id(ref_genome)
    nucl_pileup = pileup[:, :4, :]  # keep only nucleotides, no N or gap
    assert np.all(ref_genome_idxs < 4), "reference sequence includes gaps or Ns"
    # evaluate number of consensus
    n_consensus_f = np.choose(ref_genome_idxs, nucl_pileup[0])
    n_consensus_r = np.choose(ref_genome_idxs, nucl_pileup[1])
    n_cons = np.vstack([n_consensus_f, n_consensus_r])
    n_tot = np.sum(nucl_pileup, axis=1)
    return n_cons, n_tot


def safe_division(a, b, threshold=0):
    """Safe division between two arrays, returning nan if the value of the divisor
    is below threshold"""
    res = np.zeros_like(a, dtype=float) * np.nan
    mask = np.array(b) > threshold
    res = np.divide(a, b, where=mask, out=res)
    return res


def create_statsdf(Ntrue, Ntot):
    """
    Utility function to create a stats dataframe from dictionaries with the number of events.
    These dictionary have timepoints (strings) as keys. Each element is a matrix with size (2 x L)
    where L is the length of the sequence. Each row in the dataframe corresponds to one position and
    one direction (fw, rev or tot). It contains the number of observations per timepoint ("N_t")
    and the frequency of true over total observations ("freq_t") per timepoint. This is np.nan
    if no observations are present.
    """
    assert set(Ntrue.keys()) == set(
        Ntot.keys()
    ), "the two dictionaries must have the same set of keys."
    kind_dtype = pd.CategoricalDtype(["tot", "fwd", "rev"], ordered=True)
    order = np.argsort([int(t) for t in Ntrue.keys()])
    times = np.array(list(Ntrue.keys()))[order]
    P = Ntot[times[0]].shape[1]
    # initialize forward, reverse and tot dataframes
    df_f, df_r, df_t = [
        pd.DataFrame(
            {
                "position": np.arange(P, dtype=int),
                "type": pd.Series([k] * P, dtype=kind_dtype),
            }
        )
        for k in ["fwd", "rev", "tot"]
    ]
    for t in times:
        ntrue, ntot = Ntrue[t], Ntot[t]

        # Ntot
        df_f[f"N_{t}"] = ntot[0, :]
        df_r[f"N_{t}"] = ntot[1, :]
        df_t[f"N_{t}"] = ntot.sum(axis=0)

        # frequency
        df_f[f"freq_{t}"] = safe_division(ntrue[0, :], ntot[0, :])
        df_r[f"freq_{t}"] = safe_division(ntrue[1, :], ntot[1, :])
        df_t[f"freq_{t}"] = safe_division(ntrue.sum(axis=0), ntot.sum(axis=0))

    # concatenate and return
    df = pd.concat([df_t, df_f, df_r], ignore_index=True)
    df.sort_values("position", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


class StatsTable:
    """
    Class that contains a dataframe with frequency and number of observations for
    a given statistic (e.g. consensus frequency, gap density...), subdivided by
    position on the genome and forward, reverse, or total number of reads.
    """

    def __init__(self, stat, times, df):
        self.stat = stat
        self.times = times
        self.df = df

    @classmethod
    def from_counts(cls, stat_name, Ntrue, Ntot):
        """
        Creates a StatsTable given the name of the statistic and two dictionaries.
        These have timepoints as keys, and items are (2xL) matrices containing number
        of total and true observations.
        """
        stat = stat_name
        times = sorted([int(k) for k in Ntrue.keys()])
        df = create_statsdf(Ntrue, Ntot)
        return cls(stat=stat, times=times, df=df)

    def save(self, fld):
        """
        Saves the StatsTable as a pandas dataframe in the form of a gzipped picke file.
        """
        filename = pth.Path(fld) / f"stat_table_{self.stat}.pkl.gz"
        self.df.to_pickle(filename, compression="gzip")

    @classmethod
    def load(cls, filename):
        """
        Loads a StatsTable from a gzipped pickle file.
        """
        df = pd.read_pickle(filename, compression="gzip")
        stat = re.search(str(filename), r"stat_table_([^/]*).pkl.gz$").groups()[0]
        times = np.sort(
            [int(re.search(c, r"time_(\d+)$").groups()[0]) for c in list(df.columns)]
        )
        return cls(stat=stat, times=times, df=df)
