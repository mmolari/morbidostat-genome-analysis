# %%
import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

try:
    from plot_utils import *
except:
    from .plot_utils import *


# %%


def plot_1_freqhist(st, delta_freq):
    """
    Ax 1: distribution of frequencies for each time point
    Ax 2: distribution of frequency differences between final and initial time
    """
    colors = color_dict(st.times)

    cons_freqs = {t: st.freq(t, kind="tot") for t in st.times}

    fig, axs = plt.subplots(1, 2, figsize=(9, 4))

    ax = axs[0]
    bins = np.linspace(0, 1, 100)
    dict_histograms(cons_freqs, ax, colors, bins=bins)
    ax.set_yscale("log")
    ax.legend(loc="upper left")
    ax.set_xlabel("consensus frequency")
    ax.set_ylabel("n. sites")

    ax = axs[1]
    bins = np.linspace(-1, 1, 200)
    ax.hist(delta_freq, bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.set_xlabel("delta consensus frequency (final time - initial time)")
    ax.set_ylabel("n. sites")
    return fig, ax


def plot_2_pval_2dhist(pval_f, pval_r, p_thr):
    """
    Plot the joint distribution of binomial p-values for the forward and
    reverse reads, together with the marginal distributions. The red line
    represents the acceptance threshold below which positions are discarded.
    """

    fig, axs = plt.subplot_mosaic(
        """
        B.
        AC
        """,
        figsize=(8, 8),
        gridspec_kw={
            "height_ratios": [0.4, 1.0],
            "width_ratios": [1.0, 0.4],
        },
    )

    bins = np.linspace(0, 1, 50)
    norm = mpl.colors.LogNorm(vmin=0.1)

    ax = axs["A"]
    ax.hist2d(pval_f, pval_r, bins=bins, norm=norm, cmap="Blues")
    x = np.linspace(p_thr, 1, 50)
    y = p_thr / x
    ax.plot(x, y, color="red")
    ax.set_xlabel("p-value forward")
    ax.set_ylabel("p-value reverse")

    ax = axs["B"]
    ax.hist(pval_f, bins=bins)
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_ylabel("n. sites")

    ax = axs["C"]
    ax.hist(pval_r, bins=bins, orientation="horizontal")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel("n. sites")

    return fig, ax


def plot_3_pval_deltafreq(pval_tot, rank, Ne, N_thr, p_thr):
    """2D histogram of the delta-frequency vs the total p-value, separated in positions
    with more or less than `N_thr` counts."""

    mask = Ne["tot"] >= N_thr

    bins = [np.linspace(-1, 1, 50), np.linspace(0, 1, 25)]
    norm = mpl.colors.LogNorm(vmin=0.1)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for i, m in enumerate([mask, ~mask]):
        ax = axs[i]
        s = ax.hist2d(
            rank[m],
            pval_tot[m],
            bins=bins,
            norm=norm,
            cmap="Blues",
        )
        ax.axhline(p_thr, c="gray", ls="--")
        ax.set_xlabel("consensus frequency difference (final-initial)")
        ax.set_ylabel("p-value product (p-fwd*p-rev)")
        plt.colorbar(s[3], ax=ax, label="N. sites")
    axs[0].set_title(f"N $\geq$ {N_thr}")
    axs[1].set_title(f"N < {N_thr}")

    return fig, axs


def plot_4_trajectories(st, S_pos, S_rank):
    """
    For each selected position, plots the frequency trajectory.
    """

    Nkept = len(S_pos)
    Nx = 3
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)
    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    # plot one of the three trajectory lines
    def plot_single_traj(
        ax, f, n, threshold, color, marker, fillstyle="full", linestyle="-"
    ):
        if not linestyle in ["", None]:
            ax.plot(f, color=color, marker=None, linestyle=linestyle, alpha=0.3)
        kwargs = {"color": color, "fillstyle": fillstyle, "marker": marker}
        for i, fi in enumerate(f):
            if n[i] >= threshold:
                ax.plot([i], [fi], **kwargs)
            else:
                ax.plot([i], [fi], alpha=0.4, **kwargs)

    # single trajectory plot
    def plot_trajectory(ax, F, N, times, threshold=5):
        plot_single_traj(ax, F["tot"], N["tot"], threshold, "k", "o", fillstyle="none")
        plot_single_traj(ax, F["rev"], N["rev"], threshold, "C1", "+", linestyle="")
        plot_single_traj(ax, F["fwd"], N["fwd"], threshold, "C0", "x", linestyle="")

        ax.set_xticks(range(len(times)))
        ax.set_xticklabels(times)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # plot trajectories
    for ntr, p in enumerate(S_pos):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        F, N = st.traj(p)
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        plot_trajectory(ax, F, N, st.times)
        ax.set_title(f"pos = {p + 1}, $\Delta f$ = {S_rank[ntr]:.2}")

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

    fig.supylabel("consensus frequency")

    return fig, axs


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

    # data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_01")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print
    # fig_path = pth.Path("../figures/2022-05-11_RT-Tol-Res/vial_01")

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    vprint(f"preparing frequency plots for vial {vial}")

    # load consensus frequencies
    vprint("load data")
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)

    # capture values of frequencies and number of reads per site
    tb, te = st.times[0], st.times[-1]
    kinds = ["tot", "fwd", "rev"]
    Nb, Ne = [{k: st.N(t, kind=k) for k in kinds} for t in [tb, te]]
    Fb, Fe = [{k: st.freq(t, kind=k) for k in kinds} for t in [tb, te]]

    # %%
    # ~~~~~~~~~~~~ PLOT 1 ~~~~~~~~~~~~
    # Hisstogram of consensus frequencies and consensy frequency differences
    fig, ax = plot_1_freqhist(st, delta_freq=Fe["tot"] - Fb["tot"])
    plt.tight_layout()
    savefig("consensus_freq_distributions.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ FILTER ~~~~~~~~~~~~
    # rank trajectories by delta-frequency
    # evaluate binomial p-values to see if forward and reverse consensus frequencies
    # are compatible with the average frequency

    # binomial
    pval_f = binomial_pvals(Fe["fwd"], Ne["fwd"], Fe["tot"])
    pval_r = binomial_pvals(Fe["rev"], Ne["rev"], Fe["tot"])

    p_thr = 0.05
    N_thr = 10
    pval_tot = pval_f * pval_r
    rank = Fe["tot"] - Fb["tot"]
    keep = pval_tot >= p_thr
    keep &= Ne["tot"] >= N_thr
    Nkeep = 100

    # select top positions and have their ranking
    S_pos, S_rank = select_top_positions(
        ranking=rank, Nmax=Nkeep, mask=keep, inverse=True
    )

    # %%
    # ~~~~~~~~~~~~ PLOT 2 ~~~~~~~~~~~~
    # p-values forward-reverse 2d histogram

    fig, axs = plot_2_pval_2dhist(pval_f, pval_r, p_thr)
    plt.tight_layout()
    savefig("consensus_pvalue_distribution.pdf")
    show()

    # %%

    # ~~~~~~~~~~~~ PLOT 3 ~~~~~~~~~~~~
    # p-value vs delta frequency for N above or below threshold
    fig, axs = plot_3_pval_deltafreq(pval_tot, rank, Ne, N_thr, p_thr)
    plt.tight_layout()
    savefig("consensus_pvalue_vs_freqdiff.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ PLOT 4 ~~~~~~~~~~~~
    # frequency trajectories
    fig, axs = plot_4_trajectories(st, S_pos, S_rank)
    plt.tight_layout()
    savefig("consensus_min_finalfreq_trajs.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ EXPORT CSV ~~~~~~~~~~~~
    # export the selected positions as a csv dataframe with frequencies and
    # number of observations, as well as pvalues and rankings

    # add columns to the dataframe: pvalues and delta frequency
    pval_dict = {}
    rank_dict = {}
    for p, r in zip(S_pos, S_rank):
        rank_dict[p] = r
        pval_dict[(p, "fwd")] = pval_f[p]
        pval_dict[(p, "rev")] = pval_r[p]
        pval_dict[(p, "tot")] = pval_f[p] * pval_r[p]

    # select relevant positions and add columns
    mask = st.df["position"].isin(S_pos)
    sdf = st.df[mask].copy()
    sdf["pval"] = sdf.apply(lambda x: pval_dict[(x.position, x.type)], axis=1)
    sdf["rank"] = sdf.apply(lambda x: rank_dict[x.position], axis=1)

    # reorder and save
    sdf.sort_values(
        ["rank", "position", "type"], ascending=[False, True, True], inplace=True
    )
    sdf.to_csv(fig_path / "consensus_selected_positions.csv", index=False)

# %%
