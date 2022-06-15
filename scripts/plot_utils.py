import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle as pkl
import gzip
import re
import argparse

from collections import defaultdict
from Bio import SeqIO
from scipy.stats import binomtest


from extract_stats_utils import StatsTable


# Script utils
def argparser():
    """Initializes and returns an argument parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vial_fld", help="input folder containing the data for one vial"
    )
    parser.add_argument("--fig_fld", help="output folder where figure are saved")
    parser.add_argument(
        "--show", help="if specified displays the plots.", action="store_true"
    )
    parser.add_argument(
        "--verbose", help="verbose logging on stdout.", action="store_true"
    )
    return parser


# ~~~~~~~~~~~~~~~~ PROCESSING UTILS ~~~~~~~~~~~~~~~~


def average_every_step(arr, step):
    """Return an average of the array every `step` positions. Also
    returns the position corresponding to each average."""
    L = len(arr)
    N = L // step
    x = np.arange(N) * step + step / 2
    avg = arr[: N * step].reshape(-1, step).mean(axis=1)
    return x, avg


def binomial_pvals(freqs, Ns, mus):
    """Give a set of frequencies, number of observations and means, it evaluates
    the p-value for the binomial test of the frequencies relative to the fixed
    means."""
    Ks = np.rint(freqs * Ns).astype(int)
    m = np.isnan(freqs)
    pvals = np.zeros_like(freqs)
    lut = {}  # lookup table
    for i in range(len(freqs)):
        if m[i]:
            pvals[i] = np.nan
            continue
        args = (Ks[i], Ns[i], mus[i])
        if args in lut:
            pvals[i] = lut[args]
        else:
            pvals[i] = binomtest(*args).pvalue
            lut[args] = pvals[i]
    return pvals


def select_top_positions(ranking, Nmax, mask=None, inverse=False):
    """
    Select top Nmax positions in the ranking, optionally masking out values
    if `mask` argument is passed.
    """
    ranking = -np.array(ranking) if inverse else np.array(ranking)
    if mask is not None:
        ranking[~mask] = np.nan
    pos = np.arange(len(ranking))
    m = np.isnan(ranking)
    ranking, pos = ranking[~m], pos[~m]
    order = np.argsort(ranking)
    n_take = min(Nmax, len(ranking))
    idxs = order[-n_take:][::-1]
    return pos[idxs], ranking[idxs]


# ~~~~~~~~~~~~~~~~ PLOT UTILS ~~~~~~~~~~~~~~~~


def color_dict(timepoints):
    """Dictionary that assigns a color to every timepoint."""
    cmap = plt.get_cmap("rainbow")
    T = len(timepoints)
    colors = {tp: cmap(n / (T - 1)) for n, tp in enumerate(timepoints)}
    return colors


def dict_histograms(items_dict, ax, colors, plotmeans=False, **kwargs):
    """plot the histogram of a dictionary of items"""
    means = {tp: arr.mean() for tp, arr in items_dict.items()}
    for tp, arr in items_dict.items():
        ax.hist(arr, label=tp, histtype="step", color=colors[tp], **kwargs)
        if plotmeans:
            ax.axvline(means[tp], ls=":", color=colors[tp])


def plot_traj_markers(ax, f, n, thr, color, marker, fillstyle="full"):
    """
    Utility function to plot the markers of a single trajectory, The transparency
    of the marker depends on whether the number of counts is below threshold.
    """
    kwargs = {"color": color, "fillstyle": fillstyle, "marker": marker}
    for i, fi in enumerate(f):
        if n[i] >= thr:
            ax.plot([i], [fi], **kwargs)
        else:
            ax.plot([i], [fi], alpha=0.4, **kwargs)


def plot_trajectory(ax, F, N, times, thr=5, legend=False):
    """
    Given an ax, a dictionary of frequencies and number of observations, plots the
    evolution of frequencies over time. Black markers and line represent the total
    frequency, while blue and orange marker represent the forward and reverse
    frequencies respectively. Markers are transparent if the number of reads is below
    a threshold.
    """
    ax.plot(F["tot"], color="k", marker="o", alpha=0.3)
    plot_traj_markers(ax, F["tot"], N["tot"], thr, "k", "o", fillstyle="none")
    plot_traj_markers(ax, F["rev"], N["rev"], thr, "C1", "+")
    plot_traj_markers(ax, F["fwd"], N["fwd"], thr, "C0", "x")

    if legend:
        Lin = mpl.lines.Line2D
        handles = [
            Lin([0], [0], color="k", marker="o", ls="", label="tot", fillstyle="none"),
            Lin([0], [0], color="C0", marker="x", ls="", label="fwd"),
            Lin([0], [0], color="C1", marker="+", ls="", label="rev"),
        ]
        ax.legend(handles=handles, prop={"size": 6})

    ax.set_xticks(range(len(times)))
    ax.set_xticklabels(times)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ~~~~~~~~~~~~~~~~ MAIN FIGURES ~~~~~~~~~~~~~~~~


def fig_freqhist(stats_table, delta_freq):
    """
    Creates a figure with two plots. The first one is the distribution of
    frequencies over different timepoints (extracted from a StatsTable object).
    The second is the distribution of delta-frequencies. These are the difference
    between the final and initial frequency. Returns the figure and axes.
    """

    colors = color_dict(stats_table.times)
    cons_freqs = {t: stats_table.freq(t, kind="tot") for t in stats_table.times}

    fig, axs = plt.subplots(1, 2, figsize=(9, 4))

    ax = axs[0]
    bins = np.linspace(0, 1, 100)
    dict_histograms(cons_freqs, ax, colors, bins=bins)
    ax.set_yscale("log")
    ax.legend(loc="best")
    ax.set_xlabel("frequency")
    ax.set_ylabel("n. sites")

    ax = axs[1]
    bins = np.linspace(-1, 1, 200)
    ax.hist(delta_freq, bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.set_xlabel("delta frequency (final time - initial time)")
    ax.set_ylabel("n. sites")

    return fig, axs


def fig_trajectories(st, pos, rank, threshold=5):
    """
    For each selected position, plots the frequency trajectory.
    """

    Nkept = len(pos)
    Nx = 3
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)
    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    # plot trajectories
    leg_pos = min([2, Nkept])
    for ntr, p in enumerate(pos):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        F, N = st.traj(p)
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        plot_trajectory(ax, F, N, st.times, thr=threshold, legend=ntr == leg_pos)
        ax.set_title(f"pos = {p + 1}, $\Delta f$ = {rank[ntr]:.2}")

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

    return fig, axs


# def load_insertions(data_path):
#     """Given the path to the vial_X folder of interest, this function
#     loads all the corresponding insertions dictionaries. Returns them
#     in a nested dictionary whose keys are timepoint labels."""
#     # list of "time_*" folders
#     timepoint_fld = sorted(list(data_path.glob("time_*")))
#     # define regex to extract timepoint label
#     m = re.compile("/time_(.+)$")
#     # build and return dictionary of pileup matrices
#     insertions = {}
#     for tpf in timepoint_fld:
#         time_id = m.search(str(tpf)).groups()[0]
#         tp_file = tpf / "pileup" / "insertions.pkl.gz"
#         with gzip.open(tp_file, "r") as f:
#             insertions[time_id] = pkl.load(f)
#     return insertions


# def coverage(pileup):
#     """Returns the sum of nucleotide reads (no gap or N) per position"""
#     return pileup.sum(axis=0)[:4, :].sum(axis=0)


# def gap_frequency(pileup):
#     """returns the frequency of gaps (n. gaps divided by total number of reads, including N)"""
#     tot_pileup = pileup.sum(axis=0)
#     return tot_pileup[4] / tot_pileup.sum(axis=0)


# def ref_genome_to_nucl_id(ref_genome):
#     """return an array of indices (0 to 5), one per nucleotide of a
#     DNA sequence, according to the order (A,C,G,T,-,N)"""
#     nuc_alpha = ["A", "C", "G", "T", "-", "N"]
#     nuc_idx = {n: i for i, n in enumerate(nuc_alpha)}
#     nseq = np.array(ref_genome)
#     return np.array([nuc_idx[x] for x in nseq])


# def consensus_frequency(pileup, ref_genome_idxs, count_threshold=5):
#     """Returns the consensus frequencies (n. times a nucleotide different from
#     the one on the genome is observed). Gaps are discarded"""
#     tot_pileup = pileup.sum(axis=0)[:4, :]  # keep only nucleotides, no N or gap
#     assert np.all(ref_genome_idxs < 4), "reference sequence includes gaps or Ns"
#     n_consensus = np.choose(ref_genome_idxs, tot_pileup)
#     consensus_freq = np.zeros_like(n_consensus, dtype=float)
#     position_tot = tot_pileup.sum(axis=0)
#     consensus_freq = safe_division(n_consensus, position_tot, count_threshold)
#     return consensus_freq


# def consensus_frequency_from_counts(n_cons, n_tot):
#     """Returns the consensus frequencies (n. times a nucleotide different from
#     the one on the genome is observed) given the number of consensus and total
#     counts, in the form of [fwd/rev, pos] matrices."""
#     nc = n_cons.sum(axis=0)
#     nt = n_tot.sum(axis=0)
#     freq = np.zeros_like(nc, dtype=float) * np.nan
#     freq = np.divide(nc, nt, where=nt > 0, out=freq)
#     return freq


# def consensus_count(pileup, ref_genome_idxs):
#     """Returns the number of consensus nucleotides (n. times a nucleotide is
#     the one on the reference genome) and the total number of nucleotides. Gaps are discarded"""
#     nucl_pileup = pileup[:, :4, :]  # keep only nucleotides, no N or gap
#     assert np.all(ref_genome_idxs < 4), "reference sequence includes gaps or Ns"
#     # evaluate number of consensus
#     n_consensus_f = np.choose(ref_genome_idxs, nucl_pileup[0])
#     n_consensus_r = np.choose(ref_genome_idxs, nucl_pileup[1])
#     n_cons = np.vstack([n_consensus_f, n_consensus_r])
#     n_tot = np.sum(nucl_pileup, axis=1)
#     return n_cons, n_tot


# def gap_count(pileup):
#     """Given a pileup returns two matrices, one for the number of gaps and
#     one for the number of non-gap reads. The first index of the matrix refers to
#     forward or backward reads."""
#     ngaps = pileup[:, 4, :]
#     nnucl = pileup.sum(axis=1) - ngaps
#     return ngaps, nnucl


# def gap_frequency(ngaps, nnucl):
#     """Given the two matrices produced by the `gap_count` function, it evaluates
#     The gap frequency."""
#     ng = ngaps.sum(axis=0)
#     nn = nnucl.sum(axis=0)
#     return safe_division(ng, nn + ng)


# def safe_division(a, b, threshold=0):
#     """Safe division between two arrays, returning nan if the value of the divisor
#     is below threshold"""
#     res = np.zeros_like(a, dtype=float) * np.nan
#     mask = np.array(b) > threshold
#     res = np.divide(a, b, where=mask, out=res)
#     return res


# def extract_trajectories(n_true, n_tot, idxs, threshold):
#     """Extract frequency trajectories for different positions on the genome.
#     See the documentation of `extract_trajectory`"""
#     times = sorted(list(n_tot.keys()))
#     trajs = {}
#     for idx in idxs:
#         trajs[idx] = extract_trajectory(n_true, n_tot, idx, times, threshold)
#     return trajs


# def extract_trajectory(n_true, n_tot, idx, times, threshold):
#     """Given two dictionaries (n_true, n_tot) with the number of times an outcome is
#     observed and the total number of times (organized as time -> [fwd/rev , position])
#     it creates a trajectory containing the frequency of outcome on forward, reverse and total
#     observation. This is specific for the given position (idx argument).
#     The value of `times` specifies the order in which timepoints should appear. Points
#     with less than `threshold` observations are not displayed."""
#     lists = [n_true, n_tot]
#     nf, Nf = [np.array([l[t][0, idx] for t in times]) for l in lists]
#     nr, Nr = [np.array([l[t][1, idx] for t in times]) for l in lists]
#     tr = {}
#     tr["freq_f"] = safe_division(nf, Nf, threshold)
#     tr["freq_r"] = safe_division(nr, Nr, threshold)
#     tr["freq_t"] = safe_division(nf + nr, Nf + Nr, threshold)
#     return tr


# def plot_stepwise_average(arr, step, ax, **kwargs):
#     """plot the stepwise average of an array"""
#     x, avg = average_every_step(arr, step)
#     ax.plot(x, avg, **kwargs)


# # Processing Insertions


# def insertion_len(ins):
#     """Given an insertion dictionary {position -> {sequence -> [fwd, rev]}}
#     It returns the count of insertions for different lengths, in the form of
#     a dictionary {length -> n. insertions}"""
#     Ls = defaultdict(int)
#     for pos, I in ins.items():
#         for seq, N in I.items():
#             Ls[len(seq)] += np.sum(N)
#     return Ls


# def insertion_density(insertions, step, L):
#     """Given an insertion dictionary {position -> {seq -> [n. fwd, n.rev]}},
#     a step size, and a total genome length, it returns two arrays of insertion densities.
#     These arrays have a size L / step, and each contains the density of insertions in an
#     interval of length `step` basepairs.
#     The first array contains the density of insertions measured as number of insertions per position.
#     The second one as number of nucleotides inserted per position, and also includes
#     information on insertion length. These arrays are matrices with a number of rows
#     equal to L / step, and two columns. The latter correspond to forward and reverse reads."""
#     n_bins = int(np.ceil(L / step))
#     bp_density, n_density = [np.zeros((n_bins, 2), dtype=float) for _ in range(2)]
#     for pos, I in insertions.items():
#         idx = pos // step
#         for seq, n_ins in I.items():
#             n_density[idx, :] += n_ins
#             bp_density[idx, :] += len(seq) * n_ins

#     bp_density /= step
#     n_density /= step

#     L_last_chunk = L % step
#     bp_density[-1, :] *= step / L_last_chunk
#     n_density[-1, :] *= step / L_last_chunk

#     return n_density, bp_density


# def pick_insertion_threshold_length(ins_Ls, n_selected):
#     """Select a threshold length that includes the top `n_selected` longest insertions"""
#     tf = sorted(ins_Ls.keys())[-1]
#     Lf = ins_Ls[tf]
#     ins_lengths_list = sorted(Lf.keys())[::-1]
#     tot = 0
#     for l in ins_lengths_list:
#         threshold = l
#         tot += Lf[l]
#         if tot >= n_selected:
#             break
#     return threshold
