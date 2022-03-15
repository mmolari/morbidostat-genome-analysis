# %%
import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

from collections import defaultdict

try:
    from pileupplots_utils import *
except:
    from .pileupplots_utils import *


def cumulative_histograms(frequencies, ax, colors, plotmeans, **kwargs):
    """plot the cumulative histogram of frequencies"""
    means = {tp: arr.mean() for tp, arr in frequencies.items()}
    for tp, arr in frequencies.items():
        ax.hist(arr, label=tp, histtype="step", color=colors[tp], **kwargs)
        if plotmeans:
            ax.axvline(means[tp], ls=":", color=colors[tp])


def consensus_freq_signal(pileups, cons_freqs, freq_threshold, ref_genome_idxs):
    """
    Extract statistics for relevant trajectories, whose consensus frequency drops below
    the threshold at any time-point.

    Args:
        - pileups: dictionary time -> pileup matrix
        - cons_freqs: dictionary time -> cons. freq. array (one item per position)
        - freq_threshold (float): threshold consensus frequency
        - ref_genome_idx: numpy int array, containing the index (0 to 3)
            of each reference nucleotide in the pileup matrix
    Returns:
        - positions: array of relevant positions of the genome
        - freq_trajs: matrix [N. sites x timepoints] matrix with the consensus frequency of
            each relevant site for each timepoint.
        - cons_stats: dictionary containing 4 entrties. Each entry is a matrix
            of integers. Entrties are the number of consensus reads (fwd and rev)
            anf the total number of reads (fwd and rev)
    """
    # find a list of relevant positions
    positions = []
    for tp, freq in cons_freqs.items():
        x = np.arange(len(freq))
        mask = freq < freq_threshold
        positions.append(x[mask])
    positions = np.unique(np.concatenate(positions))
    # create matrix with frequency trajectories
    freq_trajs = []
    for tp, freq in cons_freqs.items():
        freq_trajs.append(freq[positions])
    freq_trajs = np.vstack(freq_trajs).T
    # create split in forward and reverse dissense from consensus
    cons_stats = defaultdict(list)
    refgen_subidxs = ref_genome_idxs[positions]
    for tp, pp in pileups.items():
        relevant_pileup = pp[:, :, positions]
        fwd_pileup = relevant_pileup[0, :4, :]
        rev_pileup = relevant_pileup[1, :4, :]
        cons_stats["ncons_fwd"].append(np.choose(refgen_subidxs, fwd_pileup))
        cons_stats["ncons_rev"].append(np.choose(refgen_subidxs, rev_pileup))
        cons_stats["ntot_fwd"].append(fwd_pileup.sum(axis=0))
        cons_stats["ntot_rev"].append(rev_pileup.sum(axis=0))

    for k in cons_stats:
        cons_stats[k] = np.vstack(cons_stats[k])

    return positions, freq_trajs, cons_stats


# %%

# path for a single vial
# input_data_path = "../results/2022-02-08_RT_test/vial_02/"

if __name__ == "__main__":

    parser = argparser()
    args = parser.parse_args()
    data_path = pth.Path(args.vial_fld)
    fig_path = pth.Path(args.fig_fld)
    fig_path.mkdir(exist_ok=True)

    # override show and save functions
    show = show_function(args.show)
    savefig = savefig_function(fig_path)

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing frequency plots for vial {vial}")

    # load pileups and reference genome
    pileups = load_pileups(data_path)
    ref_genome = load_reference_genome(data_path)

    # evaluate consensus frequencies
    ref_genome_idxs = ref_genome_to_nucl_id(ref_genome)
    cons_freqs = {
        k: consensus_frequency(pp, ref_genome_idxs) for k, pp in pileups.items()
    }

    # evaluate positions with interesting consensus frequency deviations
    freq_thr = 0.5
    end_loop = False
    while True:
        print(f"evaluating positions with cons. frequency variation, f<{freq_thr}")
        f_positions, f_traj, f_stats = consensus_freq_signal(
            pileups, cons_freqs, freq_thr, ref_genome_idxs
        )
        print(f"n. trajectories retained = {len(f_positions)}")
        if len(f_positions) <= 50:
            break

        if freq_thr <= 0.05:
            break

        freq_thr /= 1.25

    # assign colors to timepoints
    colors = color_dict(pileups)

    # %%

    # consensus frequency plot
    fig, axs = plt.subplot_mosaic("ABBB", figsize=(15, 4))
    axs["B"].set_title(f"vial {vial}")

    # frequency distribution
    ax = axs["A"]
    bins = np.linspace(0, 1, 25)
    cumulative_histograms(cons_freqs, ax, colors, plotmeans=False, bins=bins)
    ax.axvline(freq_thr, ls=":", color="k")
    ax.set_yscale("log")
    ax.set_xlabel("consensus frequency")
    ax.set_ylabel("n. sites")
    ax.legend(title="time", loc="upper left")

    # frequency vs position
    ax = axs["B"]
    for tp, freq in cons_freqs.items():
        L = freq.size
        x = np.arange(L)
        mask = freq < freq_thr
        ax.plot(x, freq, ".", label=tp, color="gray", alpha=0.05, rasterized=True)
        ax.plot(x[mask], freq[mask], ".", label=tp, color=colors[tp])
    ax.set_xlabel("position on the genome (bp)")
    ax.set_ylabel(f"consensus frequency (color if f<{freq_thr})")
    ax.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    savefig("consensus_freq_vs_position.pdf")
    show()

    # %%
    # plot frequency trajectories

    # number of trajectories and timepoints
    Ntr, T = f_traj.shape

    if Ntr > 50:
        print("WARNING: only first 50 traj plotted")
        Ntr = 50

    fig, axs = plt.subplots(Ntr, 1, figsize=(6, Ntr * 1.7), sharex=True)

    x = np.arange(T)
    t_labels = list(pileups.keys())
    for ntr in range(Ntr):
        ax = axs[ntr]
        ax.plot(x, f_traj[ntr], "o-", color="k", label="tot")
        fc, ft = f_stats["ncons_fwd"][:, ntr], f_stats["ntot_fwd"][:, ntr]
        fwd_freq = np.zeros_like(ft, dtype=float)
        fwd_freq[ft == 0] = np.nan
        fwd_freq = np.divide(fc, ft, where=ft > 0, out=fwd_freq)
        rc, rt = f_stats["ncons_rev"][:, ntr], f_stats["ntot_rev"][:, ntr]
        rev_freq = np.zeros_like(rt, dtype=float)
        rev_freq[rt == 0] = np.nan
        rev_freq = np.divide(rc, rt, where=rt > 0, out=rev_freq)
        ax.plot(x, fwd_freq, "--", marker="x", label="fwd")
        ax.plot(x, rev_freq, "--", marker="x", label="rev")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"position = {f_positions[ntr]:d}")

    ax = axs[-1]
    ax.set_xticks(x)
    ax.set_xticklabels(t_labels)
    ax.set_xlabel("timepoints")

    fig.supylabel("consensus frequency")
    plt.tight_layout()
    savefig("consensus_freq_trajectories.pdf")
    show()
