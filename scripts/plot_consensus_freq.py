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

    # override print, show, save functions
    vprint = print if args.verbose else lambda x: None
    show = plt.show if args.show else plt.close
    savefig = lambda name: plt.savefig(fig_path / name)

    # %%

    data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_02")
    savefig = lambda x: None
    show = plt.show
    vprint = print

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    vprint(f"preparing frequency plots for vial {vial}")

    vprint("load data")
    # evaluate consensus frequencies
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)

    # capture values

    tb, te = st.times[0], st.times[-1]
    kinds = ["tot", "fwd", "rev"]
    Nb, Nfb, Nrb = [st.N(tb, kind=k) for k in kinds]
    Ne, Nfe, Nre = [st.N(te, kind=k) for k in kinds]
    Fb, Ffb, Frb = [st.freq(tb, kind=k) for k in kinds]
    Fe, Ffe, Fre = [st.freq(te, kind=k) for k in kinds]

    # %%
    # evaluate quantitites of interest and plots
    pval_f = binomial_pvals(Ffe, Nfe, Fe)
    pval_r = binomial_pvals(Fre, Nre, Fe)

    rank = Fe - Fb
    keep = (pval_f > 0.05) & (pval_r > 0.05)
    keep &= Fb >= 0.8
    Nkeep = 100

    pos, r = select_top_positions(ranking=rank, Nmax=Nkeep, mask=keep, inverse=True)

    # %%
    # CONSENSUS PLOT 1)
    # Histogram of frequencies and of difference of frequencies between initial and final points
    colors = color_dict(st.times)

    # plot frequency trajectories

    cons_freqs = {t: st.freq(t, kind="tot") for t in st.times}
    # plot frequency of gaps
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))

    ax = axs[0]
    bins = np.linspace(0, 1, 100)
    dict_histograms(cons_freqs, ax, colors, bins=bins)
    ax.set_yscale("log")
    ax.legend(loc="upper left")
    ax.set_xlabel("consensus frequency")
    ax.set_ylabel("n. sites")

    ax = axs[1]
    bins = np.linspace(-1, 1, 200)
    ax.hist(Fe - Fb, bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.set_xlabel("delta consensus frequency (final time - initial time)")
    ax.set_ylabel("n. sites")

    plt.tight_layout()
    savefig("consensus_freq_histograms.pdf")
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

# TODO:
# - csv file with positions, delta frequency and p value
# - plot two histograms
# - plot trajectories
# - pvalue plots with thresholds
# - document plots
# - more agressive filtering of pvalues (product? Both more than 50%?)
# - scatter plot pvalue and delta frequency
