# %%
import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

try:
    from utils.plot_utils import *
except:
    from .utils.plot_utils import *


def plot_n_ins_genome(st, step):
    """
    Plots the frequency of insertions, for different timepoints and for fwd/rev reads,
    cumulative over a window of size `step` bps.
    """

    Ts = np.sort(st.times)
    NT = len(Ts)

    fig, axs = plt.subplots(NT, 1, figsize=(12, NT * 4), sharex=True)

    for nt, t in enumerate(Ts):
        ax = axs[nt]

        pos = np.sort(st.df.position.unique())
        bins = np.arange(0, pos.max() + 2 * step, step)

        for k in ["fwd", "rev"]:
            g = 1 if k == "fwd" else -1
            F = st.freq(t, k)
            F = np.nan_to_num(F, nan=0.0)
            ax.hist(
                pos,
                weights=F * g,
                bins=bins,
                histtype="step",
                label=k,
                rasterized=True,
            )
        ax.legend()
        ax.set_xlim(bins[0] - step * 5, bins[-1] + step * 5)
        ax.set_xlabel("genome position (bp)")
        ax.set_ylabel(f"n. gaps per {step} bp")
        ax.grid(alpha=0.3)
        ax.set_title(f"t = {t}")
        ytl = [int(y) for y in ax.get_yticks()]
        ax.set_yticks(ytl)
        ax.set_yticklabels([str(y).lstrip("-") for y in ytl])

        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1e6))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1e5))
        ax.grid(alpha=0.2, which="major")
        ax.grid(alpha=0.1, which="minor")

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

    # data_path = pth.Path("../results/2022-amoxicilin/vial_7")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print
    # fig_path = pth.Path("../figures/2022-05-11_RT-Tol-Res/vial_01")

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    vprint(f"preparing gap frequency plots for vial {vial}")

    # load gap frequencies
    vprint("load data")
    st_path = data_path / "stats" / "stats_table_gap_freq.pkl.gz"
    st = StatsTable.load(st_path)

    # capture values of frequencies and number of reads per site
    tb, te = st.times[0], st.times[-1]
    kinds = ["tot", "fwd", "rev"]
    Nb, Ne = [{k: st.N(t, kind=k) for k in kinds} for t in [tb, te]]
    Fb, Fe = [{k: st.freq(t, kind=k) for k in kinds} for t in [tb, te]]

    # %%
    # ~~~~~~~~~~~~ PLOT 1 ~~~~~~~~~~~~
    # Histogram of gap frequencies and consensy frequency differences
    vprint(f"preparing plot 1: gap frequency distribution")

    fig, axs = fig_freqhist(st, delta_freq=Fe["tot"] - Fb["tot"])

    axs[0].set_xlabel("gap frequency")
    axs[1].set_xlabel("delta gap frequency (final time - initial time)")

    plt.tight_layout()
    savefig("gap_freq_distributions.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ PLOT 2 ~~~~~~~~~~~~
    # gap frequency histogram over time (fwd/rev)
    fig, axs = plot_n_ins_genome(st, step=100)
    plt.tight_layout()
    savefig("gap_vs_genome.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ FILTER ~~~~~~~~~~~~
    # rank trajectories by delta-frequency
    # evaluate binomial p-values to see if forward and reverse gap frequencies
    # are compatible with the average frequency
    vprint(f"selecting positions with high delta frequency")

    p_thr = 0.05
    N_thr = 10

    rank, keep = st.rank_trajs(p_threshold=p_thr, n_threshold=N_thr)

    Nkeep = 100

    # select top positions and have their ranking
    S_pos, S_rank = select_top_positions(
        ranking=rank, Nmax=Nkeep, mask=keep, inverse=False
    )

    # %%
    # ~~~~~~~~~~~~ PLOT 4 ~~~~~~~~~~~~
    # frequency trajectories
    vprint(f"preparing plot 4: frequency trajectories of selected sites")

    fig, axs = fig_trajectories(st, S_pos, S_rank)
    fig.supylabel("gap frequency")
    fig.supxlabel("timepoint")
    plt.tight_layout()
    savefig("gap_freq_trajs.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ EXPORT CSV ~~~~~~~~~~~~
    # export the selected positions as a csv dataframe with frequencies and
    # number of observations, as well as pvalues and rankings
    vprint(f"export csv file with selected positions")

    # add columns to the dataframe: pvalues and delta frequency
    rank_dict = {}
    for p, r in zip(S_pos, S_rank):
        rank_dict[p] = r

    # select relevant positions and add columns
    mask = st.df["position"].isin(S_pos)
    sdf = st.df[mask].copy()
    sdf["rank"] = sdf.apply(lambda x: rank_dict[x.position], axis=1)

    # make position 1-based
    sdf["position"] += 1

    # reorder and save
    sdf.sort_values(
        ["rank", "position", "type"], ascending=[False, True, True], inplace=True
    )
    sdf.to_csv(fig_path / "gap_selected_positions.csv", index=False)

# %%
