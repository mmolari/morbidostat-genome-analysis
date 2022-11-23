# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

import pandas as pd

from matplotlib.ticker import MultipleLocator

try:
    from utils.plot_utils import *
    from utils.extract_stats_utils import safe_division
except:
    from .utils.plot_utils import *
    from .utils.extract_stats_utils import safe_division


# %%


def build_clip_dataframes(clips, clip_seqs):
    """Function that builds a set of clip dataframes, given the
    list of clips and the StatsTable object containing the consensus
    frequency.
    The return value is a dictionary with times as keys and dataframes as
    values. The dataframes contain 4 x 3 columns, corresponding to number
    of clips (C), total number of read start/ends at position (N), frequency of
    clips (F) and average length of soft clips (L) for forward, reverse
    and total reads (f,r,t).
    """

    dfs = {}
    for t, clps in clips.items():

        pos = np.array(sorted(clps.keys()))
        df = {}

        # n. clips
        C = np.vstack([clps[p] for p in pos])
        df["Cf"], df["Cr"] = C[:, 0], C[:, 1]
        df["Ct"] = df["Cf"] + df["Cr"]

        # number of reads
        df["Nf"], df["Nr"] = C[:, 2], C[:, 3]
        df["Nt"] = df["Nf"] + df["Nr"]

        # frequency of clips
        df["Ff"] = safe_division(df["Cf"], df["Nf"])
        df["Fr"] = safe_division(df["Cr"], df["Nr"])
        df["Ft"] = safe_division(df["Ct"], df["Nt"])

        # build dataframe
        df = pd.DataFrame(df, index=pos)

        # average length of clipped part
        df["Lf"], df["Lr"], df["Lt"] = np.nan, np.nan, np.nan

        cs = clip_seqs[t]
        pf, pr, pt = [], [], []
        Lf, Lr, Lt = [], [], []
        for p, v in cs.items():
            lf = [len(s) for s in v[0]]
            lr = [len(s) for s in v[1]]
            if len(lf) > 0:
                Lf.append(np.mean(lf))
                pf.append(p)
            if len(lr) > 0:
                Lr.append(np.mean(lr))
                pr.append(p)
            Lt.append(np.mean(lf + lr))
            pt.append(p)

        for lab, p, val in zip(["Lf", "Lr", "Lt"], [pf, pr, pt], [Lf, Lr, Lt]):
            df.loc[p, lab] = val

        dfs[t] = df

    return dfs


def plot_histograms(dfs, Ts):
    """Plots the histograms of number of clips, clip frequencies and
    average clipped length (Ct, Ft, Lt) for all the different timepoints."""

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    colors = color_dict(Ts)
    ax = axs[0]
    dx = {k: df["Ct"] for k, df in dfs.items()}
    M = max([df.max() for df in dx.values()])
    bins = np.arange(1, M + 2) - 0.5
    dict_histograms(dx, ax, colors, bins=bins)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("n. clips")
    ax.set_ylabel("n. sites")
    ax.legend()

    ax = axs[1]
    dx = {k: df["Ft"] for k, df in dfs.items()}
    dict_histograms(dx, ax, colors, bins=np.linspace(0, 1, 100))
    ax.set_yscale("log")
    ax.set_xlabel("clip frequency")
    ax.set_ylabel("n. sites")

    ax = axs[2]
    dx = {k: df["Lt"] for k, df in dfs.items()}
    M = max([df.max() for df in dx.values()])
    m = max([df.min() for df in dx.values()])
    bins = np.logspace(np.log10(m) - 0.1, np.log10(M) + 0.1, 50)
    dict_histograms(dx, ax, colors, bins=bins)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("average soft-clipped length (bp)")
    ax.set_ylabel("n. sites")

    return fig, axs


def plot_joint_distr(df):
    """Plots the joint distribution of forward and reverse values for
    the number of clips, clip frequency and clipped length
    (C, F, L) corresponding to the last timepoint."""

    fig, axs = plt.subplots(1, 3, figsize=(14, 4))

    ax = axs[0]
    f, r = df["Cf"], df["Cr"]
    bins = np.arange(max([f.max(), r.max()]) + 2) - 0.5
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("n. clips forward")
    ax.set_ylabel("n. clips reverse")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    ax = axs[1]
    f, r = df["Ff"], df["Fr"]
    bins = np.linspace(0, 1, 25)
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("freq clips forward")
    ax.set_ylabel("freq clips reverse")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    ax = axs[2]
    f, r = df["Lf"], df["Lr"]
    M = max([f.max(), r.max()])
    m = max([f.min(), r.min()])
    bins = np.logspace(np.log10(m) - 0.1, np.log10(M) + 0.1, 25)
    norm = mpl.colors.LogNorm()
    m = ax.hist2d(f, r, bins=bins, norm=norm)
    ax.set_xlabel("avg. clip length forward")
    ax.set_ylabel("avg. clip length reverse")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.colorbar(m[3], ax=ax, label="n. sites")

    return fig, axs


def plot_n_clips_genome(dfs, step=5000):
    """For each timepoint plots the histogram of number of clips per genome position.
    The number is the sum over intervals of length `step`. It is reported separately
    for forward and reverse reads."""

    Ts = sorted(dfs.keys())
    NT = len(Ts)

    fig, axs = plt.subplots(NT, 1, figsize=(12, NT * 4), sharex=True)

    for nt, t in enumerate(Ts):
        ax = axs[nt]
        df = dfs[t]
        sf, sr = df["Cf"], df["Cr"]
        M = max([sf.index.max(), sr.index.max()])
        bins = np.arange(0, M + step, step)
        for s, l in zip([sf, sr], ["fwd", "rev"]):
            g = 1 if l == "fwd" else -1
            ax.hist(s.index, weights=s.values * g, bins=bins, histtype="step", label=l)
        ax.legend()
        ax.set_xlim(bins[0] - step * 5, bins[-1] + step * 5)
        ax.set_xlabel("genome position (bp)")
        ax.set_ylabel(f"n. clips per {step//1000} kbp")
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


def plot_trajectories(sel_df, Ts, threshold=5, Nx=3):
    """Plot the trajectories of n. of clips for the selected trajectories.
    Black lines represent the total n. clips, and blue and orange
    markers represent respectively forward and reverse n. clips."""

    Nkept = len(sel_df)
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 2.0)

    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    kinds = ["tot", "fwd", "rev"]

    # plot trajectories
    Tf = Ts[-1]
    leg_pos = min([2, Nkept])
    for ntr, it in enumerate(sel_df.iterrows()):
        p, row = it
        axidx = np.unravel_index(ntr, (Ny, Nx))
        C = {k: [row[f"C{k[0]}_{t}"] for t in Ts] for k in kinds}
        N = {k: [row[f"N{k[0]}_{t}"] for t in Ts] for k in kinds}
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        plot_trajectory(ax, C, N, Ts, thr=threshold, legend=ntr == leg_pos)
        Ctots = [row[f"Ct_{t}"] for t in Ts]
        idx_max = np.argmax(Ctots)
        Cmax = int(Ctots[idx_max])
        Tmax = Ts[idx_max]
        l = row[f"Lt_{Tmax}"]
        ax.set_title(f"{p+1}, Cmax={Cmax}, " + r"$\langle L \rangle$" + f"={int(l)}")
        ax.grid(alpha=0.3)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_ylim(bottom=0)

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

    fig.supylabel("n. clips")
    fig.supxlabel("timepoint")

    return fig, axs


def select_clip_positions(dfs, Nkeep):
    """Select clip positions with maximum delta-count over time."""

    # list of all positions over all timepoints
    clip_positions = [df.index.to_numpy() for df in dfs.values()]
    clip_positions = np.sort(np.unique(np.concatenate(clip_positions)))
    pos_idx = {p: i for i, p in enumerate(clip_positions)}

    # list of times
    Ts = list(sorted(dfs.keys()))

    # create empty count container
    Npos = len(clip_positions)
    Ntimes = len(Ts)
    Nclips = np.zeros((Ntimes, Npos), dtype=int)

    # fill with number of clip counts per position
    for n, t in enumerate(Ts):
        df = dfs[t]
        C = df["Ct"]
        idxs = [pos_idx[p] for p in df.index]
        Nclips[n, idxs] = C

    # evaluate delta = max - min
    Nclips = Nclips.T
    delta = np.max(Nclips, axis=1) - np.min(Nclips, axis=1)

    # select positions with maximum delta
    order = np.argsort(delta)[::-1]
    idx_keep = order[:Nkeep]
    selected_pos = np.sort(clip_positions[idx_keep])

    return selected_pos


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

    # debug purpose
    # data_path = pth.Path("../results/2022-amoxicilin/vial_7")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print
    # fig_path = pth.Path("../figures/2022-05-11_RT-Tol-Res/vial_04")

    # %%

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing clip plots for vial {vial}")

    # load clips and consensus frequency
    clips, clip_seqs = load_clips(data_path)
    Ts = sorted(clips.keys())
    Tf = Ts[-1]  # final time

    # %%
    # build clip dataframes
    dfs = build_clip_dataframes(clips, clip_seqs)

    # %%
    # plot histogram of N. clips, clip frequency, clip length
    fig, axs = plot_histograms(dfs, Ts)
    plt.tight_layout()
    savefig("clips_distributions.pdf")
    show()

    # %%
    # joint distribution for the last timepoint
    df = dfs[Tf]
    fig, axs = plot_joint_distr(df)
    plt.tight_layout()
    savefig("clips_joint_dist_final_time.pdf")
    show()

    # %%
    # histogram of number of clips per different timepoints
    fig, axs = plot_n_clips_genome(dfs)
    plt.tight_layout()
    savefig("clips_vs_genome.pdf")
    show()

    # %%
    # select positions with maximum delta-number of clips (WARNING: not normalized)
    Nkeep = 51  # number of max selected positions
    selected_pos = select_clip_positions(dfs, Nkeep)

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
    savefig("clips_trajectories.pdf")
    show()

    # %%
    # format dataframe and export to csv:
    # separate each row of sel_df in three rows (tot, fwd, rev) and save
    # the entries as F, N, C, L for each timepoint.

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
    exp_df.to_csv(fig_path / "clips_selected_positions.csv", index=False)

# %%
