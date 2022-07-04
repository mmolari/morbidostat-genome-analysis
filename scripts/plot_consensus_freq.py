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
    vprint(f"preparing consensus frequency plots for vial {vial}")

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
    # Histogram of consensus frequencies and consensy frequency differences
    vprint(f"preparing plot 1: consensus frequency distribution")

    fig, axs = fig_freqhist(st, delta_freq=Fe["tot"] - Fb["tot"])

    axs[0].set_xlabel("consensus frequency")
    axs[1].set_xlabel("delta consensus frequency (final time - initial time)")

    plt.tight_layout()
    savefig("consensus_freq_distributions.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ FILTER ~~~~~~~~~~~~
    # rank trajectories by delta-frequency
    # evaluate binomial p-values to see if forward and reverse consensus frequencies
    # are compatible with the average frequency
    vprint(f"selecting positions with high delta frequency")

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
    vprint(f"preparing plot 2: p-value joint distribution")

    fig, axs = plot_2_pval_2dhist(pval_f, pval_r, p_thr)
    plt.tight_layout()
    savefig("consensus_pvalue_distribution.pdf")
    show()

    # %%

    # ~~~~~~~~~~~~ PLOT 3 ~~~~~~~~~~~~
    # p-value vs delta frequency for N above or below threshold
    vprint(f"preparing plot 3: p-value vs delta frequency distribution")

    fig, axs = plot_3_pval_deltafreq(pval_tot, rank, Ne, N_thr, p_thr)
    plt.tight_layout()
    savefig("consensus_pvalue_vs_freqdiff.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ PLOT 4 ~~~~~~~~~~~~
    # frequency trajectories
    vprint(f"preparing plot 4: frequency trajectories of selected sites")

    fig, axs = fig_trajectories(st, S_pos, S_rank)
    fig.supylabel("consensus frequency")
    fig.supxlabel("timepoint")
    plt.tight_layout()
    savefig("consensus_freq_trajs.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ EXPORT CSV ~~~~~~~~~~~~
    # export the selected positions as a csv dataframe with frequencies and
    # number of observations, as well as pvalues and rankings
    vprint(f"export csv file with selected positions")

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

    # make position 1-based
    sdf["position"] += 1

    # reorder and save
    sdf.sort_values(
        ["rank", "position", "type"], ascending=[False, True, True], inplace=True
    )
    sdf.to_csv(fig_path / "consensus_selected_positions.csv", index=False)

# %%
