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


def plot_n_gaps_vs_genome(st, step):
    """
    Plots the cumulative frequency of gaps for different timepoints and for
    fwd/rev reads, over a window of size `step` bps.
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


class glob_interval:
    """Class that describes an interval of positions, used to
    group trajectories together."""

    def __init__(self, p):
        self.beg = p
        self.end = p
        self.len = 1

    def __lt__(self, other):
        """sort based on beginning position"""
        return self.beg < other.beg

    def extend(self, q):
        """Extend the interval in one direction"""
        if q == self.beg - 1:
            self.beg -= 1
            self.len += 1
            return True
        elif q == self.end + 1:
            self.end += 1
            self.len += 1
            return True
        else:
            return False

    def attempt_merge(self, other):
        """attempt to merge two indervals if beg-end are only 1 apart.
        Modifies self and returns True in case of success."""
        if self.end == other.beg - 1:
            self.end = other.end
            return True
        elif self.beg == other.end + 1:
            self.beg = other.beg
            return True
        return False

    def __str__(self):
        return f"interval [{self.beg},{self.end}]"

    def to_list(self):
        """Returns a list of interval positions"""
        return list(range(self.beg, self.end + 1))


def select_and_glob_best_trajs(rank, keep, Nkeep):
    """Select Nkeep trajectory groups, returning a list of positions for each group.
    Trajectories are ordered based on their rank, and only trajectores with
    keep=True are considered.
    Trajectories are progressively added to the list of groups. If a position falls
    on the edge of a group, it is added to the group instead of forming a new group.
    Groups can be merged together if they have matching edges. This stops once
    Nkeep groups are formed and no more trajectories can be added to the groups
    without creating a new one."""

    # order trajectories by rank
    order = np.argsort(rank)[::-1]

    # list of intervals to populate
    intervals = []

    for o in order:

        # skeep trajectories with keep=False
        if not keep[o]:
            continue

        # try to extend an existing interval
        extend = False
        for i in intervals:
            extend |= i.extend(o)
            if extend:
                break

        if not extend:
            # if failure, create a new singleton group with the position
            if len(intervals) == Nkeep:
                break
            intervals.append(glob_interval(o))
        else:
            # if success, check if two adjacent intervals can be merged
            intervals = list(sorted(intervals))
            for j in range(len(intervals) - 1):
                merged = intervals[j].attempt_merge(intervals[j + 1])
                if merged:
                    intervals.pop(j + 1)
                    break
    # return a list of positions for each interval
    return [i.to_list() for i in intervals]


def fig_glob_trajectories(st, G_pos, G_rank, threshold=5):
    """
    For each selected position group, plots the frequency trajectory.
    """

    # # reorder based on rank
    # Avg_rank = [np.mean(gr) for gr in G_rank]
    # order = np.argsort(Avg_rank)[::-1]
    # G_pos = [G_pos[o] for o in order]
    # G_rank = [G_rank[o] for o in order]

    Nkept = len(G_pos)
    Nx = 3
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)
    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    # plot trajectories
    leg_pos = min([2, Nkept])
    for ntr, pos in enumerate(G_pos):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        for p in pos:
            F, N = st.traj(p)
            plot_trajectory(ax, F, N, st.times, thr=threshold, legend=ntr == leg_pos)
        avg_rank = np.mean(G_rank[ntr])
        if len(pos) == 1:
            ax.set_title(f"{pos[0] + 1}, $\Delta f$ = {avg_rank:.2}")
        else:
            ax.set_title(
                rf"[{pos[0] + 1},{pos[-1] + 1}], $\langle \Delta f \rangle$ = {avg_rank:.2}"
            )

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

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
    fig, axs = plot_n_gaps_vs_genome(st, step=100)
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

    Nkeep = 99

    # select top positions (glob adjacent positions together)
    G_pos = select_and_glob_best_trajs(rank, keep, Nkeep)
    G_rank = [[rank[p] for p in i] for i in G_pos]

    # %%
    # ~~~~~~~~~~~~ PLOT 4 ~~~~~~~~~~~~
    # frequency trajectories
    vprint(f"preparing plot 4: frequency trajectories of selected sites")

    fig, axs = fig_glob_trajectories(st, G_pos, G_rank)
    # fig.supylabel("gap frequency")
    # fig.supxlabel("timepoint")
    plt.tight_layout()
    savefig("gap_freq_trajs.pdf")
    show()

    # %%
    # ~~~~~~~~~~~~ EXPORT CSV ~~~~~~~~~~~~
    # export the selected positions as a csv dataframe with frequencies and
    # number of observations, as well as pvalues and rankings
    vprint(f"export csv file with selected positions")

    # add columns to the dataframe: pvalues and delta frequency
    S_pos = np.concatenate(G_pos)
    S_rank = np.concatenate(G_rank)
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
