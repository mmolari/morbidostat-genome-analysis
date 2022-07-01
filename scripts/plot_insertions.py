# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

import pandas as pd

try:
    from plot_utils import *
    from extract_stats_utils import safe_division
except:
    from .plot_utils import *
    from .extract_stats_utils import safe_division


# %%


def n_ins(x):
    return np.vstack(list(x.values())).sum(axis=0)


def L_tot(x):
    return np.vstack([I * len(seq) for seq, I in x.items()]).sum(axis=0)


def build_insertion_dataframes(insertions, stats_table):
    """Function that builds a set of insertion dataframes, given the
    list of insertions and the StatsTable object containing the consensus
    frequency.
    The return value is a dictionary with times as keys and dataframe
    """

    Ts = stats_table.times

    dfs = {}
    for t in Ts:

        ins = insertions[str(t)]
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
    data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_04")
    savefig = lambda x: None
    show = plt.show
    vprint = print
    fig_path = pth.Path("../figures/2022-05-11_RT-Tol-Res/vial_04")

    # %%

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing insertion plots for vial {vial}")

    # load insertions and consensus frequency
    insertions = load_insertions(data_path)
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)

    # %%

    # build insertion dataframe
    dfs = build_insertion_dataframes(insertions, st)

    # %%
    # plot histogram of N. insertions, insertion frequency, insertion length

    colors = color_dict(st.times)
    # mask_f = lambda df: (df["Ff"] > 0) & (df["Fr"] > 0) & (df["It"] > 2)

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

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

    plt.tight_layout()
    savefig("insertions_distributions.pdf")
    show()

    # %%
    # joint distribution

    df = dfs[st.times.max()]

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

    plt.tight_layout()
    savefig("insertions_joint_dist_final_time.pdf")
    show()

    # %%
    # histogram of number of insertions
    df = dfs[st.times.max()]

    fig, axs = plt.subplots(2, 1, figsize=(12, 8))

    ax = axs[0]
    sf, sr = df["If"], df["Ir"]
    M = max([sf.index.max(), sr.index.max()])
    step = 5000
    bins = np.arange(0, M + step, step)
    ax.hist(sf.index, weights=sf.values, bins=bins, histtype="step", label="fwd")
    ax.hist(sr.index, weights=sr.values, bins=bins, histtype="step", label="rev")
    ax.legend()
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel(f"n. insertions per {step//1000} kbp")

    ax = axs[1]
    s = df["Lt"]
    M = s.index.max()
    step = 5000
    bins = np.arange(0, M + step, step)
    ax.hist(s.index, weights=s.values * df["Ft"], bins=bins, histtype="step", color="k")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel(f"avg. n. inserted bp per {step//1000} kbp")

    plt.tight_layout()
    savefig("insertions_vs_genome.pdf")
    show()

    # %%
    # select relevant positions

    df = dfs[st.times.max()]

    selected_pos = []
    NposI, NposF, NposL = 40, 40, 20

    # 1) select positions with highest total number of insertions
    # exclude positions with less than one insertion in both directions
    mask = (df["If"] > 0) & (df["Ir"] > 0)
    selected_pos += list(df[mask]["It"].sort_values(ascending=False)[:NposI].index)

    # 2) select position with high insertions frequency in both fwd and reverse (max fwd + rev)
    selected_pos += list(df[mask]["Ft"].sort_values(ascending=False)[:NposF].index)

    # 3) select positions with long insertions
    selected_pos += list(df[mask]["Lt"].sort_values(ascending=False)[:NposL].index)

    selected_pos = np.sort(np.unique(selected_pos))

    # %%
    # create dataframe
    Tf = st.times.max()
    sel_df = dfs[Tf].loc[selected_pos].add_suffix(f"_{Tf}")
    for t in st.times[::-1]:
        if t == Tf:
            continue
        sel_df = sel_df.join(dfs[t].add_suffix(f"_{t}"), how="left")
    sel_df = sel_df.fillna(0)
    tp = lambda k: float if k.startswith("F") else int
    sel_df = sel_df.astype({k: tp(k) for k in sel_df.columns})

    # %%
    # draw trajectories

    threshold = 4
    Nkept = len(selected_pos)
    Nx = 3
    Ny = int(np.ceil(Nkept / Nx))
    figsize = (Nx * 3, Ny * 1.5)
    fig, axs = plt.subplots(Ny, Nx, figsize=figsize, sharex=True, sharey=True)

    kinds = ["tot", "fwd", "rev"]

    # plot trajectories
    Ts = np.sort(st.times)
    Tf = Ts[-1]
    leg_pos = min([2, Nkept])
    for ntr, p in enumerate(selected_pos):
        axidx = np.unravel_index(ntr, (Ny, Nx))
        row = sel_df.loc[p]
        F = {k: [row[f"F{k[0]}_{t}"] for t in Ts] for k in kinds}
        N = {k: [row[f"N{k[0]}_{t}"] for t in Ts] for k in kinds}
        axidx = axidx[1] if Ny == 1 else axidx
        ax = axs[axidx]
        plot_trajectory(ax, F, N, st.times, thr=threshold, legend=ntr == leg_pos)
        i = int(row[f"It_{Tf}"])
        l = row[f"Lt_{Tf}"]
        ax.set_title(f"pos={p+1}, I={i}, " + r"$\langle L \rangle$" + f"= {l:.2}")

    # remove extra axes
    for i in range(Nkept, Nx * Ny):
        axidx = np.unravel_index(i, (Ny, Nx))
        axidx = axidx[1] if Ny == 1 else axidx
        axs[axidx].remove()

    plt.tight_layout()
    savefig("insertions_trajectories.pdf")
    show()

    # %%

    # export dataframe to csv

    # make positions 1-based
    sel_df.index += 1

    kinds = ["tot", "fwd", "rev"]
    exp_df = []
    for p, row in sel_df.iterrows():
        for k in kinds:
            it = {"position": p, "kind": k}
            for l, v in row.iteritems():
                if l[1] == k[0]:
                    ln = l[0] + l[-2:]
                    it[ln] = v
            exp_df.append(it)
    exp_df = pd.DataFrame(exp_df)
    exp_df.to_csv(fig_path / "insertions_selected_positions.csv", index=False)

# %%
