# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

try:
    from plot_utils import *
except:
    from .plot_utils import *

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
    print(f"preparing gap plots for vial {vial}")

    # load pileups and reference genome
    pileups = load_pileups(data_path)

    # get gap counts
    ngaps, nnucl, ntot = {}, {}, {}
    for k, pp in pileups.items():
        ngaps[k], nnucl[k] = gap_count(pp)
        ntot[k] = ngaps[k] + nnucl[k]

    # evaluate frequency of gaps
    gap_freq = {k: gap_frequency(ngaps[k], nnucl[k]) for k in ngaps}

    # evaluate delta frequencies
    times = sorted(list(gap_freq.keys()))
    ti, tf = times[0], times[-1]
    delta_freq = gap_freq[tf] - gap_freq[ti]

    # number of positions to keep
    N_keep = 50

    # find positions with highest gap final frequency
    top_idx, top_vals = find_top_N(gap_freq[tf], N_keep)
    freq_thr = top_vals[0]

    # find positions with highest delta gap frequency
    topdf_idx, top_dfvals = find_top_N(delta_freq, N_keep)
    dfreq_thr = top_dfvals[0]

    # extract frequency trajectories
    n_threshold = 5
    traj_f = extract_trajectories(ngaps, ntot, top_idx, n_threshold)
    traj_df = extract_trajectories(ngaps, ntot, topdf_idx, n_threshold)

    colors = color_dict(gap_freq)
    # %%

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
    cumulative_histograms(gap_freq, ax, colors, plotmeans=False, bins=bins)
    ax.set_yscale("log")
    ax.axvline(freq_thr, ls="--", color="k")
    ax.legend(loc="upper right")
    ax.set_xlabel("gap frequency")
    ax.set_ylabel("n. sites")

    ax = axs["B"]
    for k, freq in gap_freq.items():
        x = np.arange(len(freq))
        ax.plot(x[top_idx], freq[top_idx], ".", color=colors[k])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel("gap frequency")

    ax.set_xlim(0, x.max())

    ax = axs["C"]
    bins = np.linspace(-1, 1, 200)
    ax.hist(delta_freq, bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.axvline(dfreq_thr, ls="--", color="gray")
    ax.set_xlabel("delta gap frequency (final time - initial time)")
    ax.set_ylabel("n. sites")

    ax = axs["D"]
    x = np.arange(len(freq))
    ax.plot(x[topdf_idx], delta_freq[topdf_idx], "k.")
    ax.set_xlim(0, x.max())
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel("delta gap frequency (tf - ti)")

    plt.tight_layout()
    savefig("gap_freq_vs_position.pdf")
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
    fig.supylabel("gap frequency")
    plt.tight_layout()
    savefig("gap_max_finalfreq_trajs.pdf")
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
    fig.supylabel("gap frequency")
    plt.tight_layout()
    savefig("gap_max_deltafreq_trajs.pdf")
    show()
# %%
