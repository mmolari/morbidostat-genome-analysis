# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
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

    # data_path = pth.Path("../results/2022-05-11_RT-Tol-Res/vial_02")
    # savefig = lambda x: None
    # show = plt.show
    # vprint = print

    # extract vial
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    vprint(f"preparing coverage plots for vial {vial}")

    vprint(f"loading data and extracting coverages")
    # load data
    st_path = data_path / "stats" / "stats_table_reference_freq.pkl.gz"
    st = StatsTable.load(st_path)
    coverages = {t: st.N(t, kind="tot") for t in st.times}

    # assign colors to timepoints
    colors = color_dict(st.times)

    # %%
    # COVERAGE PLOT 1)
    # left part: cumulative histogram of coverage per site
    # Right part: average coverage over the sequence

    vprint("plotting coverage distribution")
    # figure setup
    fig, axs = plt.subplot_mosaic("ABBBB", figsize=(25, 4))

    # histogram of coverages
    ax = axs["A"]
    maxcov = max([max(cov) for cov in coverages.values()]) + 2
    bins = np.arange(maxcov)
    dict_histograms(coverages, ax, colors, plotmeans=True, bins=bins, cumulative=True)
    ax.set_xlim(right=maxcov - 1)
    ax.legend(loc="lower right", title="time")
    ax.set_xlabel("coverage per site")
    ax.set_ylabel("n. sites (cumulative)")

    # coverage along the genome (mean of every kbp)
    ax = axs["B"]
    ax.set_title(f"vial {vial}")
    step = 1000
    for tp, cov in coverages.items():
        cov_m = cov / cov.mean()
        x, subcov = average_every_step(cov_m, step)
        ax.plot(x, subcov, color=colors[tp], alpha=0.8)

    # guidelines (mean and 2 std)
    cov_mat = np.vstack(list(coverages.values()))
    cov_rescaled = (cov_mat.T / cov_mat.mean(axis=1)).T
    cov_mean = cov_rescaled.mean(axis=0)
    cov_std = cov_rescaled.std(axis=0)
    m_step = 50000
    xm, cov_mean = average_every_step(cov_mean, m_step)
    xs, cov_std = average_every_step(cov_std, m_step)
    ax.plot(
        xm,
        cov_mean,
        "--",
        color="gray",
        label=f"mean ({m_step//1000} kbp) " + r"$\pm 2 \sigma$",
    )
    ax.plot(xm, cov_mean - 2 * cov_std, ":", color="gray")
    ax.plot(xm, cov_mean + 2 * cov_std, ":", color="gray")

    # axes setup
    ax.legend(loc="lower right")
    ax.set_xlabel("position on the genome (bp)")
    ax.set_ylabel(f"coverage ({step//1000} kbp mean) / mean coverage")
    ax.set_xlim(0, xm.max())
    plt.tight_layout()
    savefig("coverage_distribution.pdf")
    show()

    # %%
    # COVERAGE PLOT 2)
    # forward and reverse coverage boxplot for different timepoints

    vprint("plotting forward/reverse coverage boxplot")

    # create dataframe with forward and reverse coverage
    dfs = []
    for k in ["fwd", "rev"]:
        for t in st.times:
            df = pd.DataFrame(pd.Series(st.N(t, kind=k), name="coverage"))
            df["orientation"] = k
            df["time"] = t
            dfs.append(df)
    dfs = pd.concat(dfs)

    # draw the boxplot
    fig, ax = plt.subplots(1, 1, figsize=(len(st.times) + 1, 4))
    sns.boxplot(
        data=dfs, x="time", y="coverage", hue="orientation", showfliers=False, ax=ax
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig("coverage_fwd_rev.pdf")
    show()

# %%
