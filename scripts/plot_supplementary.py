import argparse

import pandas as pd
import matplotlib.pyplot as plt

from itertools import combinations
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that plots a matrix + histogram describing positions of suppl. alignments.
    """
    )
    parser.add_argument("--df", type=str, help="non_primary.csv dataframe")
    parser.add_argument("--pdf", type=str, help="name of output pdf file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # load dataframe
    df = pd.read_csv(args.df)

    # select supplementary and primary reads
    sdf = df[df.suppl]
    pdf = df[(~df.sec) & (~df.suppl)].set_index("read", verify_integrity=True)

    # find stetches to link
    links = []

    # group by read id
    for read_id, rdf in sdf.groupby("read"):

        # find and add corresponding primary read
        primary = pdf.loc[[read_id]]
        rdf = pd.concat([rdf, primary.reset_index()])

        # select paris of reads with the same read-id
        N = rdf.shape[0]
        for i, j in combinations(range(N), 2):
            a, b = rdf.iloc[i], rdf.iloc[j]

            # check whether the reads are non-overlapping
            a_bef_b = a.qe <= b.qs
            b_bef_a = a.qs >= b.qe

            # if non-overlapping, define contact point and add to links
            l1, l2 = None, None
            if a_bef_b:
                l1 = a.re if a.fwd else a.rs
                l2 = b.rs if b.fwd else b.re
            elif b_bef_a:
                l1 = a.rs if a.fwd else a.re
                l2 = b.re if b.fwd else b.rs
            if l1 is not None:
                same_str = a.fwd == b.fwd
                links.append((l1, l2, same_str))

    # Perform the plot
    fig, axs = plt.subplots(
        2, 1, figsize=(9, 11), gridspec_kw={"height_ratios": [1, 0.1]}, sharex=True
    )
    ax = axs[0]

    # scatter plot
    X, Y, C = [[l[i] for l in links] for i in range(3)]
    C = ["k" if c else "r" for c in C]
    ax.scatter(X, Y, alpha=0.03, color=C, rasterized=True)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label=k, markerfacecolor=c, alpha=0.3)
        for k, c in zip(["fwd", "rev"], ["k", "r"])
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    ax.axis("equal")
    ax.set_xlabel("primary read location (bp)")
    ax.set_ylabel("suppl. read location (bp)")

    # diagonal and ax limits
    X = sdf.re.max()
    ax.plot([0, X], [0, X], ls=":", color="gray")
    ax.set_xlim(-X * 0.02, X * 1.02)
    ax.set_ylim(-X * 0.02, X * 1.02)

    for k in [ax.xaxis, ax.yaxis]:
        k.set_major_locator(MultipleLocator(1e6))
        k.set_minor_locator(MultipleLocator(1e5))

    ax.grid(alpha=0.2, which="major")
    ax.grid(alpha=0.1, which="minor")

    ax = axs[1]
    ax.hist(sdf.rs, bins=1000, color="k")

    ax.grid(alpha=0.2, which="major")
    ax.grid(alpha=0.1, which="minor")

    ax.set_xlabel("primary read location (bp)")
    ax.set_ylabel("n. suppl. reads")

    plt.tight_layout()
    plt.savefig(args.pdf)
    plt.close(fig)
