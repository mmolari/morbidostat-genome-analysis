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

    # capture values
    tb, te = st.times[0], st.times[-1]
    kinds = ["tot", "fwd", "rev"]
    Nb, Ne = [{k: st.N(t, kind=k) for k in kinds} for t in [tb, te]]
    Fb, Fe = [{k: st.freq(t, kind=k) for k in kinds} for t in [tb, te]]

    # %%
    # PLOT !) Hisstogram of consensus frequencies and consensy frequency differences
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
    ax.hist(Fe["tot"] - Fb["tot"], bins=bins, histtype="step", color="k")
    ax.set_yscale("log")
    ax.set_xlabel("delta consensus frequency (final time - initial time)")
    ax.set_ylabel("n. sites")

    plt.tight_layout()
    savefig("consensus_freq_distributions.pdf")
    show()

    # %%
    # FILTER
    # rank trajectories by delta-frequency
    # evaluate binomial p-values to see if forward and reverse consensus frequencies
    # are compatible with the average frequency

    pval_f = binomial_pvals(Fe["fwd"], Ne["fwd"], Fe["tot"])
    pval_r = binomial_pvals(Fe["rev"], Ne["rev"], Fe["tot"])

    p_thr = 0.05
    rank = Fe["tot"] - Fb["tot"]
    keep = pval_f * pval_r >= p_thr
    keep &= Ne["tot"] >= 10
    keep &= Fb["tot"] >= 0.6
    keep &= Fe["tot"] <= 0.85
    Nkeep = 100

    # select top positions and have their ranking
    S_pos, S_rank = select_top_positions(
        ranking=rank, Nmax=Nkeep, mask=keep, inverse=True
    )

    # %%
    # PLOT 2) p-values scatter
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

    ax = axs["A"]
    b = ax.hist2d(pval_f, pval_r, bins=bins, norm=mpl.colors.LogNorm(), cmap="Blues")
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

    plt.tight_layout()
    savefig("consensus_pvalue_distribution.pdf")
    show()

    # %%
    # PLOT 3) frequency trajectories

    Nkept = len(S_pos)
    Nx = 3
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)
    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

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

    def plot_trajectory(ax, F, N, times, threshold=5):
        plot_single_traj(ax, F["tot"], N["tot"], threshold, "k", "o", fillstyle="none")
        plot_single_traj(ax, F["rev"], N["rev"], threshold, "C1", "+", linestyle="")
        plot_single_traj(ax, F["fwd"], N["fwd"], threshold, "C0", "x", linestyle="")

        ax.set_xticks(range(len(times)))
        ax.set_xticklabels(times)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for ntr, p in enumerate(S_pos):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        F, N = st.traj(p)
        ax = axs[axidx]
        plot_trajectory(ax, F, N, st.times)
        ax.set_title(f"position = {p}")

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axs[axidx].remove()

    fig.supylabel("consensus frequency")
    plt.tight_layout()
    savefig("consensus_min_finalfreq_trajs.pdf")
    show()

    # %%
    # EXTRACT CSV
    # export the selected positions as a csv dataframe with frequencies and
    # number of observations, as well as pvalues and rankings

    pval_dict = {}
    rank_dict = {}
    for p, r in zip(S_pos, S_rank):
        rank_dict[p] = r
        pval_dict[(p, "fwd")] = pval_f[p]
        pval_dict[(p, "rev")] = pval_r[p]
        pval_dict[(p, "tot")] = pval_f[p] * pval_r[p]

    mask = st.df["position"].isin(S_pos)
    sdf = st.df[mask].copy()
    sdf["pval"] = sdf.apply(lambda x: pval_dict[(x.position, x.type)], axis=1)
    sdf["rank"] = sdf.apply(lambda x: rank_dict[x.position], axis=1)
    sdf.sort_values(
        ["rank", "position", "type"], ascending=[False, True, True], inplace=True
    )
    sdf.to_csv(fig_path / "consensus_selected_positions.csv", index=False)
# %%
# TODO:
# - document plots
# - scatter plot pvalue and delta frequency
