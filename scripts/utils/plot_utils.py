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


from .extract_stats_utils import StatsTable


# ~~~~~~~~~~~~~~~~ SCRIPT UTILS ~~~~~~~~~~~~~~~~


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


# ~~~~~~~~~~~~~~~~ LOADING UTILS ~~~~~~~~~~~~~~~~


def load_time_dict(data_path, filename):
    """Given the path to the vial_X folder of interest, this function
    loads a dictionary per timepoint, returning it in a nested
    dictionary whose keys are timepoint labels.
    The dictionaries are named `time_X/pileup/filename`."""
    # list of "time_*" folders
    timepoint_fld = sorted(list(data_path.glob("time_*")))
    # define regex to extract timepoint label
    m = re.compile("/time_(.+)$")
    # build and return dictionary of pileup matrices
    dicts = {}
    for tpf in timepoint_fld:
        time_id = m.search(str(tpf)).groups()[0]
        time_id = int(time_id)
        tp_file = tpf / "pileup" / filename
        with gzip.open(tp_file, "r") as f:
            dicts[time_id] = pkl.load(f)
    return dicts


load_insertions = lambda data_path: load_time_dict(data_path, "insertions.pkl.gz")


def load_clips(data_path):
    """Given the path to the vial_X folder of interest, this function
    loads the dictionaries of clip count and clip sequences for all timepoints.
    It returns them as two nested dictionaries whose keys are timepoints."""
    dicts = load_time_dict(data_path, "clips.pkl.gz")

    clip_count, clip_seqs = {}, {}

    for t, d in dicts.items():
        clip_count[t] = d["count"]
        clip_seqs[t] = d["seqs"]

    return clip_count, clip_seqs


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
    ordered_Ts = sorted(items_dict.keys())
    for tp in ordered_Ts:
        arr = items_dict[tp]
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
