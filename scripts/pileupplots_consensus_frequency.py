# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

try:
    from pileupplots_utils import *
except:
    from .pileupplots_utils import *


# %%

if __name__ == "__main__":

    parser = argparser()
    args = parser.parse_args()
    data_path = pth.Path(args.vial_fld)
    fig_path = pth.Path(args.fig_fld)
    fig_path.mkdir(exist_ok=True)

    # override show and save functions
    show = show_function(args.show)
    savefig = savefig_function(fig_path)

    # %%

    # data_path = pth.Path("../results/2022-02-08_RT_test/vial_04/")
    # savefig = lambda x: None
    # show = lambda: plt.show()

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing frequency plots for vial {vial}")

    # load pileups and reference genome
    pileups = load_pileups(data_path)
    ref_genome = load_reference_genome(data_path)

    # evaluate number of times a nucleotide is consensus and total number of nucleotides:
    ref_genome_idxs = ref_genome_to_nucl_id(ref_genome)
    ncons, ntot = {}, {}
    for k, pp in pileups.items():
        ncons[k], ntot[k] = consensus_count(pp, ref_genome_idxs)

    # evaluate consensus frequencies
    cons_freqs = {k: consensus_frequency_from_counts(ncons[k], ntot[k]) for k in ncons}

    # delta consensus frequency
    times = sorted(list(cons_freqs.keys()))
    ti, tf = times[0], times[-1]
    delta_freq = cons_freqs[tf] - cons_freqs[ti]

    # %%
    # number of positions to keep
    N_keep = 50

    # find positions with highest gap final frequency
    top_idx, top_vals = find_bottom_N(cons_freqs[tf], N_keep)
    freq_thr = top_vals[-1]

    # find positions with highest delta gap frequency
    topdf_idx, top_dfvals = find_bottom_N(delta_freq, N_keep)
    dfreq_thr = top_dfvals[-1]

    # extract frequency trajectories
    n_threshold = 10
    traj_f = extract_trajectories(ncons, ntot, top_idx, n_threshold)
    traj_df = extract_trajectories(ncons, ntot, topdf_idx, n_threshold)

    # assign colors to timepoints
    colors = color_dict(pileups)

    # %%
    # plot frequency trajectories

    # plot frequency of gaps
    fig, axs = plt.subplot_mosaic(
        """
        ABBB
        CDDD
        """,
        figsize=(15, 7),
    )

    ax = axs["A"]
    bins = np.linspace(0, 1, 100)
    cumulative_histograms(cons_freqs, ax, colors, plotmeans=False, bins=bins)
    ax.set_yscale("log")
    ax.axvline(freq_thr, ls="--", color="k")
    ax.legend(loc="upper right")
    ax.set_xlabel("consensus frequency")
    ax.set_ylabel("n. sites")

    ax = axs["B"]
    for k, freq in cons_freqs.items():
        x = np.arange(len(freq))
        ax.plot(x[top_idx], freq[top_idx], ".", color=colors[k])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel("consensus frequency")

    ax.set_xlim(0, x.max())

    ax = axs["C"]
    bins = np.linspace(-1, 1, 200)
    ax.hist(delta_freq, bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.axvline(dfreq_thr, ls="--", color="gray")
    ax.set_xlabel("delta consensus frequency (final time - initial time)")
    ax.set_ylabel("n. sites")

    ax = axs["D"]
    x = np.arange(len(freq))
    ax.plot(x[topdf_idx], delta_freq[topdf_idx], "k.")
    ax.set_xlim(0, x.max())
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel("delta consensus frequency (tf - ti)")

    plt.tight_layout()
    savefig("consensus_freq_vs_position.pdf")
    show()
    # %%

    Nx = 3
    Ny = int(np.ceil(N_keep / Nx))

    fig, axs = plt.subplots(
        Ny, Nx, figsize=(Nx * 3, Ny * 1.5), sharex=True, sharey=True
    )

    for ntr, pos in enumerate(traj_f):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        traj = traj_f[pos]
        ax = axs[axidx]
        plot_trajectory(ax, traj)
        ax.set_title(f"position = {pos}")

    for ax in axs[-1, :]:
        ax.set_xticks(np.arange(len(times)))
        ax.set_xticklabels(times)
    fig.supylabel("consensus frequency")
    plt.tight_layout()
    savefig("consensus_min_finalfreq_trajs.pdf")
    show()

    # %%

    fig, axs = plt.subplots(
        Ny, Nx, figsize=(Nx * 3, Ny * 1.5), sharex=True, sharey=True
    )

    for ntr, pos in enumerate(traj_df):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        traj = traj_df[pos]
        ax = axs[axidx]
        plot_trajectory(ax, traj)
        ax.set_title(f"position = {pos}")

    for ax in axs[-1, :]:
        ax.set_xticks(np.arange(len(times)))
        ax.set_xticklabels(times)
    fig.supylabel("consensus frequency")
    plt.tight_layout()
    savefig("consensus_min_deltafreq_trajs.pdf")
    show()
# %%
