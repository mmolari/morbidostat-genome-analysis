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

    Ts = stats_table.times

    dfs = {}
    for t in Ts:

        ins = insertions[str(t)]
        pos = np.array(sorted(ins.keys()))
        df = {}

        # insertions
        I = np.vstack([n_ins(ins[p]) for p in pos])
        df["If"], df["Ir"] = I[:, 0], I[:, 1]
        df["I"] = I.sum(axis=1)

        # number of reads
        df["Nf"] = stats_table.N(t, kind="fwd")[pos]
        df["Nr"] = stats_table.N(t, kind="rev")[pos]
        df["N"] = df["Nf"] + df["Nr"]

        # frequency of insertions
        df["Ff"] = safe_division(df["If"], df["Nf"])
        df["Fr"] = safe_division(df["Ir"], df["Nr"])
        df["F"] = safe_division(df["I"], df["N"])

        # average read length
        Ltot = np.vstack([L_tot(ins[p]) for p in pos])
        df["Lf"] = safe_division(Ltot[:, 0], df["If"])
        df["Lr"] = safe_division(Ltot[:, 1], df["Ir"])
        df["L"] = Ltot.sum(axis=1) / df["I"]

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

    data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_04")
    savefig = lambda x: None
    show = plt.show
    vprint = print
    fig_path = pth.Path("../figrues/2022-05-11_RT-Tol-Res/vial_04")

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
    colors = color_dict(st.times)

    for t in st.times:
        df = dfs[t]
        vals = df["I"]
        bins = np.arange(vals.max() + 2) - 0.5
        plt.hist(vals, label=f"t = {t}", histtype="step", bins=bins, color=colors[t])
    plt.yscale("log")
    plt.legend()
    plt.show()

    bins = np.linspace(0, 1, 100)
    for t in st.times:
        df = dfs[t]
        plt.hist(df["F"], label=f"t = {t}", histtype="step", bins=bins, color=colors[t])
    plt.yscale("log")
    plt.legend()
    plt.show()

# %%
