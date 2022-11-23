# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

import pandas as pd

from matplotlib.ticker import MultipleLocator
from functools import cache
from scipy.stats import fisher_exact
from collections import defaultdict


try:
    from utils.plot_utils import *
    from utils.extract_stats_utils import safe_division
except:
    from .utils.plot_utils import *
    from .utils.extract_stats_utils import safe_division


# %%


def build_insertion_dataframes(insertions, stats_table):
    """Function that builds a set of insertion dataframes, given the
    list of insertions and the StatsTable object containing the consensus
    frequency.
    The return value is a dictionary with times as keys and dataframes as
    values. The dataframes contain 4 x 3 columns, corresponding to number
    of insertions (I), total number of reads (N) frequency of insertions
    on total number of treads (F) and average length of insertions (L)
    for forward, reverse and total reads (f,r,t).
    """

    def n_ins(x):
        return np.vstack(list(x.values())).sum(axis=0)

    def L_tot(x):
        return np.vstack([I * len(seq) for seq, I in x.items()]).sum(axis=0)

    Ts = stats_table.times

    dfs = {}
    for t in Ts:

        ins = insertions[t]
        pos = np.array(sorted(ins.keys()))
        df = {}

        # insertions
        I = np.vstack([n_ins(ins[p]) for p in pos])
        df["If"], df["Ir"] = I[:, 0], I[:, 1]
        df["It"] = I.sum(axis=1)

        # number of reads
        df["Nf"] = stats_table.N(t, kind="fwd")[pos]
        df["Nr"] = stats_table.N(t, kind="rev")[pos]
        df["Nt"] = df["Nf"] + df["Nr"]

        # frequency of insertions
        df["Ff"] = safe_division(df["If"], df["Nf"])
        df["Fr"] = safe_division(df["Ir"], df["Nr"])
        df["Ft"] = safe_division(df["It"], df["Nt"])

        # average read length
        Ltot = np.vstack([L_tot(ins[p]) for p in pos])
        df["Lf"] = safe_division(Ltot[:, 0], df["If"])
        df["Lr"] = safe_division(Ltot[:, 1], df["Ir"])
        df["Lt"] = Ltot.sum(axis=1) / df["It"]

        # build dataframe
        dfs[t] = pd.DataFrame(df, index=pos)
    #     df = pd.concat(columns, axis=1).fillna(0).astype(int)

    return dfs


def plot_histograms(dfs, Ts):
    """Plots the histograms of number of insertions, insertion
    frequencies and average length (It, Ft, Lt) for all the different
    timepoints."""

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    colors = color_dict(Ts)
    ax = axs[0]
    dx = {k: df["It"] for k, df in dfs.items()}
    M = max([df.max() for df in dx.values()])
    bins = np.arange(M + 2) - 0.5
    dict_histograms(dx, ax, colors, bins=bins)
    ax.set_yscale("log")
    ax.set_xlabel("n. insertions")
    ax.set_ylabel("n. sites")
    ax.legend()

    ax = axs[1]
    dx = {k: df["Ft"] for k, df in dfs.items()}
    dict_histograms(dx, ax, colors, bins=np.linspace(0, 1, 100))
    ax.set_yscale("log")
    ax.set_xlabel("insertion frequency")
    ax.set_ylabel("n. sites")

    ax = axs[2]
    dx = {k: df["Lt"] for k, df in dfs.items()}
    M = max([df.max() for df in dx.values()])
    bins = np.logspace(0, np.log10(M), 50)
    dict_histograms(dx, ax, colors, bins=bins)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("insertion average length (bp)")
    ax.set_ylabel("n. sites")

    return fig, axs


def plot_joint_distr(df):
    """Plots the joint distribution of forward and reverse values for
    the insertion number, insertion frequency and insertion length
    (I, F, L) corresponding to the last timepoint."""

    fig, axs = plt.subplots(1, 3, figsize=(14, 4))

    ax = axs[0]
    f, r = df["If"], df["Ir"]
    bins = np.arange(max([f.max(), r.max()]) + 2) - 0.5
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("n. insertions forward")
    ax.set_ylabel("n. insertions reverse")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    ax = axs[1]
    f, r = df["Ff"], df["Fr"]
    bins = np.linspace(0, 1, 25)
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("freq insertions forward")
    ax.set_ylabel("freq insertions reverse")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    ax = axs[2]
    f, r = df["Lf"], df["Lr"]
    M = max([f.max(), r.max()])
    bins = np.logspace(0, np.log10(M), 25)
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("avg. insertion length forward")
    ax.set_ylabel("avg. insertion length reverse")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    return fig, axs


def plot_n_ins_genome(dfs, stat, step=5000):
    """
    depending on the value of `stats`, plots:
    `N`: the number of insertions in forward and reverse reads over
         a window of size `step` bps.
    `L`: the total length of insertions over all the reads on a window
        of size `step` bps.
    """

    Ts = sorted(dfs.keys())
    NT = len(Ts)

    fig, axs = plt.subplots(NT, 1, figsize=(12, NT * 4), sharex=True)

    for nt, t in enumerate(Ts):
        ax = axs[nt]
        df = dfs[t]
        if stat == "N":
            sf, sr = df["If"], df["Ir"]
        elif stat == "L":
            sf = (df["Ff"] * df["Lf"]).fillna(0)
            sr = (df["Fr"] * df["Lr"]).fillna(0)
        else:
            raise ValueError("stat must be either 'N' or 'L'")
        M = max([sf.index.max(), sr.index.max()])
        bins = np.arange(0, M + step, step)
        for s, l in zip([sf, sr], ["fwd", "rev"]):
            g = 1 if l == "fwd" else -1
            ax.hist(
                s.index,
                weights=s.values * g,
                bins=bins,
                histtype="step",
                label=l,
                rasterized=True,
            )
        ax.legend()
        ax.set_xlim(bins[0] - step * 5, bins[-1] + step * 5)
        ax.set_xlabel("genome position (bp)")
        if stat == "N":
            ax.set_ylabel(f"n. insertions per {step} bp")
        elif stat == "L":
            ax.set_ylabel(f"avg. bp inserted per {step} bp")
        ax.grid(alpha=0.3)
        ax.set_title(f"t = {t}")
        ytl = [int(y) for y in ax.get_yticks()]
        ax.set_yticks(ytl)
        ax.set_yticklabels([str(y).lstrip("-") for y in ytl])

        ax.xaxis.set_major_locator(MultipleLocator(1e6))
        ax.xaxis.set_minor_locator(MultipleLocator(1e5))
        ax.grid(alpha=0.2, which="major")
        ax.grid(alpha=0.1, which="minor")

    return fig, axs


@cache
def trust_trajectory(Tf, Nf, Tr, Nr, p_min=0.05, Nmin=10):
    """Performs an exact fisher thest on a single trajectory timepoint to
    decide whether to trust it. This function is cached.
    Returns a pair (True/False, Freq. true)"""
    if (Nf < Nmin) or (Nr < Nmin):
        return False, None
    Ff, Fr = Nf - Tf, Nr - Tr
    _, p = fisher_exact([[Tf, Ff], [Tr, Fr]])
    if p < p_min:
        return False, None
    return True, (Tf + Tr) / (Nf + Nr)


def complete_noins_trajs(trajs, Ts, pos, st, Nmin):
    """Given the trajectory frequency matrix, it completes tests for the positions with
    no insertions. These are the positions with zero frequency. For these positions the
    number of reads is extracted from the stats-table, and the function checks which sites
    have above-threshold number of reads (N >= Nmin). The rest of the frequencies are set
    to np.nan."""
    for nt, t in enumerate(Ts):
        for kind in ["fwd", "rev"]:
            mask = trajs[nt, :] == 0
            Ns = st.N(t, kind=kind)[pos]
            mask &= Ns < Nmin
            trajs[nt, mask] = np.nan
    return trajs


def select_relevant_positions(dfs, st, Nselect, Nmin, p_min):
    """Given the dictionary of dataframes, selects positions with maximal insertion
    frequency variation (max - min over timepoints).
    The variation is only considered over timepoints with more than Nmin reads in the pileup,
    and where a fisher exact test on the fwd/rev insertion frequency returns a
    p value > threshold.
    Returns a list of positions."""

    # sorted list of timepoints
    Ts = list(sorted(dfs.keys()))
    Nt = len(Ts)

    # list of all positions in which an insertion is detected at least once in a timepoint
    positions = np.sort(
        np.unique(np.concatenate([list(df.index) for df in dfs.values()]))
    )
    Np = len(positions)

    # create empty frequency trajectory matrix
    trajs = np.zeros((Nt, Np))

    # dictionary position -> index in thrajectory matrix
    pos_idxs = {p: i for i, p in enumerate(positions)}

    for nt, t in enumerate(Ts):
        print(f"processing time {t}")
        df = dfs[t]

        for pos, row in df.iterrows():

            if pos % 100000 == 0:
                print(trust_trajectory.cache_info())

            # fisher test on whether to trust the trajectory. Cached.
            Nf, Nr, If, Ir = row[["Nf", "Nr", "If", "Ir"]]
            trust, F = trust_trajectory(If, Nf, Ir, Nr, p_min, Nmin)

            pos_i = pos_idxs[pos]
            if trust:
                trajs[nt, pos_i] = F
            else:
                trajs[nt, pos_i] = np.nan

    # check time-position pairs without insertions.
    trajs = complete_noins_trajs(trajs, Ts, positions, st, Nmin)

    # rank trajectories based on max-min trusted frequencies
    trajs = trajs.T
    ranks = np.nanmax(trajs, axis=1) - np.nanmin(trajs, axis=1)

    # mask out positions with no trusted point and sort them
    mask = ~np.isnan(ranks)
    R = ranks[mask]
    P = positions[mask]
    order = np.argsort(R)[::-1]

    # select top ones
    selected_pos = P[order[:Nselect]]
    return selected_pos


def plot_trajectories(sel_df, Ts, threshold=5, Nx=3):
    """Plot the frequency trajectories (n. insertions / tot n. reads) for the selected
    trajectories. Black lines represent the total frequency, and blue and orange
    markers represent respectively forward and reverse frequencies."""

    Nkept = len(sel_df)
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)

    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    kinds = ["tot", "fwd", "rev"]

    # plot trajectories
    Ts = np.sort(st.times)
    Tf = Ts[-1]
    leg_pos = min([2, Nkept])
    for ntr, it in enumerate(sel_df.iterrows()):
        p, row = it
        axidx = np.unravel_index(ntr, (Ny, Nx))
        F = {k: [row[f"F{k[0]}_{t}"] for t in Ts] for k in kinds}
        N = {k: [row[f"N{k[0]}_{t}"] for t in Ts] for k in kinds}
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        plot_trajectory(ax, F, N, Ts, thr=threshold, legend=ntr == leg_pos)
        i = int(row[f"It_{Tf}"])
        l = row[f"Lt_{Tf}"]
        ax.set_title(f"pos={p+1}, I={i}, " + r"$\langle L \rangle$" + f"={l:.2}")

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

    fig.supylabel("insertion frequency")
    fig.supxlabel("timepoint")

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

    # # debug purpose
    # data_path = pth.Path("../results/2022-amoxicilin/vial_7")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print
    # fig_path = pth.Path("../figures/2022-amoxicilin/vial_7")

    # %%

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing insertion plots for vial {vial}")

    # load insertions and consensus frequency
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)  # consensus frequency for each site
    Ts = np.sort(st.times)  # list of times
    Tf = Ts[-1]  # final time

    # %%
    # build insertion dataframe
    dfs = build_insertion_dataframes(load_insertions(data_path), st)

    # %%
    # plot histogram of N. insertions, insertion frequency, insertion length

    fig, axs = plot_histograms(dfs, Ts)
    plt.tight_layout()
    savefig("insertions_distributions.pdf")
    show()

    # %%
    # joint distribution for the last timepoint

    df = dfs[Tf]
    fig, axs = plot_joint_distr(df)
    plt.tight_layout()
    savefig("insertions_joint_dist_final_time.pdf")
    show()

    # %%
    # histogram of number and length of insertions vs genome

    fig, axs = plot_n_ins_genome(dfs, stat="N", step=100)
    plt.tight_layout()
    savefig("insertions_vs_genome_N.pdf")
    show()

    fig, axs = plot_n_ins_genome(dfs, stat="L", step=100)
    plt.tight_layout()
    savefig("insertions_vs_genome_L.pdf")
    show()

    # %%
    # select relevant positions

    Nselect = 81
    Nmin = 8
    p_min = 0.05
    selected_pos = select_relevant_positions(dfs, st, Nselect, Nmin, p_min)

    # %%
    # create dataframe with only selected positions
    sel_df = pd.DataFrame([], index=selected_pos)
    for t in Ts[::-1]:
        sel_df = sel_df.join(dfs[t].add_suffix(f"_{t}"), how="left")
    sel_df = sel_df.fillna(0)
    tp = lambda k: float if k.startswith("F") else int
    sel_df = sel_df.astype({k: tp(k) for k in sel_df.columns})

    # %%
    # draw trajectories
    fig, axs = plot_trajectories(sel_df, Ts, threshold=5)
    plt.tight_layout()
    savefig("insertions_trajectories.pdf")
    show()

    # %%
    # format dataframe and export to csv:
    # separate each row of sel_df in three rows (tot, fwd, rev) and save
    # the entries as F, N, I, L for each timepoint.

    kinds = ["tot", "fwd", "rev"]
    exp_df = []
    for p, row in sel_df.iterrows():
        for k in kinds:
            it = {"position": p + 1, "kind": k}
            for l, v in row.iteritems():
                if l[1] == k[0]:
                    ln = l[0] + l[-2:]
                    it[ln] = v
            exp_df.append(it)
    exp_df = pd.DataFrame(exp_df)
    exp_df.to_csv(fig_path / "insertions_selected_positions.csv", index=False)

# %%
