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
    sdf.to_csv(fig_path / "consensus_selected_positions.csv", index=False)

# %%
