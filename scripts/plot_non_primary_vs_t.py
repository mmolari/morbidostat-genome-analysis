import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that plots the evolution of the histogram of secondary/supplementary
    reads over time.
    """
    )
    parser.add_argument("--dfs", nargs="+", help="list of non_primary.csv dataframe")
    parser.add_argument("--ts", type=int, nargs="+", help="list of timepoints")
    parser.add_argument("--pdf_sec", type=str, help="output pdf for secondary reads")
    parser.add_argument(
        "--pdf_suppl", type=str, help="output pdf for supplementary reads"
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # load dataframes
    dfs = [pd.read_csv(df) for df in args.dfs]
    ts = args.ts

    assert len(dfs) == len(ts)
    N = len(ts)

    # define dataframes as dictionary
    dfs = {t: df for t, df in zip(ts, dfs)}

    # define common binning
    M = max([df.re.max() for df in dfs.values()])
    bins = np.linspace(0, M, 1000)

    # ensure ordering is correct
    ts = np.sort(ts)

    # mask functions
    only_sec = lambda df: df[df.sec]
    only_suppl = lambda df: df[df.suppl]

    # two histograms, one for supplementary and one for secondary reads
    for only_f, fname, kind in zip(
        [only_sec, only_suppl],
        [args.pdf_sec, args.pdf_suppl],
        ["secondary", "supplementary"],
    ):

        fig, axs = plt.subplots(N, 1, sharex=True, figsize=(10, N * 2.3))
        for n, t in enumerate(ts):
            ax = axs[n]

            # select only secondary or supplementary reads
            df = only_f(dfs[t])

            # histogram of read starts
            ax.hist(df.rs, bins=bins, color="k")

            # set major and minor axis locator
            ax.xaxis.set_major_locator(MultipleLocator(1e6))
            ax.xaxis.set_minor_locator(MultipleLocator(1e5))

            ax.grid(alpha=0.2, which="major")
            ax.grid(alpha=0.1, which="minor")

            ax.set_title(f"t = {t}")
            ax.set_ylabel(f"n. {kind} reads")

        axs[-1].set_xlabel("reference position (bp)")

        plt.tight_layout()
        plt.savefig(fname)
        plt.close(fig)
