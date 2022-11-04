# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

import pandas as pd

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


def plot_n_ins_genome(df, step=5000):
    """Plots:
    1) the number of insertions in forward and reverse reads over
        (step) bp windows.
    2) the average number of inserted basepairs over (step) bp windows. It weights
        the insertion length by insertion frequency."""

    fig, axs = plt.subplots(2, 1, figsize=(12, 8))

    ax = axs[0]
    sf, sr = df["If"], df["Ir"]
    M = max([sf.index.max(), sr.index.max()])
    bins = np.arange(0, M + step, step)
    for s, l in zip([sf, sr], ["fwd", "rev"]):
        ax.hist(s.index, weights=s.values, bins=bins, histtype="step", label=l)
    ax.legend()
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel(f"n. insertions per {step//1000} kbp")
    ax.grid(alpha=0.3)

    ax = axs[1]
    s = df["Lt"]
    M = s.index.max()
    bins = np.arange(0, M + step, step)
    ax.hist(s.index, weights=s.values * df["Ft"], bins=bins, histtype="step", color="k")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel(f"avg. n. inserted bp per {step//1000} kbp")
    ax.grid(alpha=0.3)

    return fig, axs


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

    # debug purpose
    # data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_04")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print
    # fig_path = pth.Path("../figures/2022-05-11_RT-Tol-Res/vial_04")

    # %%

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing insertion plots for vial {vial}")

    # load insertions and consensus frequency
    insertions = load_insertions(data_path)
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)  # consensus frequency for each site
    Ts = np.sort(st.times)  # list of times
    Tf = Ts[-1]  # final time

    # %%
    # build insertion dataframe
    dfs = build_insertion_dataframes(insertions, st)

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
    # histogram of number of insertions

    df = dfs[Tf]
    fig, axs = plot_n_ins_genome(df)
    plt.tight_layout()
    savefig("insertions_vs_genome.pdf")
    show()

    # %%
    # select relevant positions

    df = dfs[Tf]

    # For most of the selected trajectories, exclude positions having no
    # insertions on forward or reverse reads
    mask = (df["If"] > 0) & (df["Ir"] > 0)

    selected_pos = []
    NposI, NposF, NposL = 40, 40, 20  # number of max selected positions (may coincide)

    # 1) select positions with highest total number of insertions
    selected_pos += list(df["It"].sort_values(ascending=False)[:3].index)
    selected_pos += list(df[mask]["It"].sort_values(ascending=False)[:NposI].index)

    # 2) select position with high insertions frequency in both fwd and reverse (max fwd + rev)
    selected_pos += list(df["Ft"].sort_values(ascending=False)[:3].index)
    selected_pos += list(df[mask]["Ft"].sort_values(ascending=False)[:NposF].index)

    # 3) select positions with long insertions
    selected_pos += list(df["Lt"].sort_values(ascending=False)[:3].index)
    selected_pos += list(df[mask]["Lt"].sort_values(ascending=False)[:NposL].index)

    selected_pos = np.sort(np.unique(selected_pos))

    # %%
    # create dataframe with only selected positions
    sel_df = dfs[Tf].loc[selected_pos].add_suffix(f"_{Tf}")
    for t in Ts[::-1]:
        if t == Tf:
            continue
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
