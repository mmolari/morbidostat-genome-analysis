import argparse

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that plots a matrix + histogram describing positions of secondary alignments.
    """
    )
    parser.add_argument("--df", type=str, help="non_primary.csv dataframe")
    parser.add_argument("--pdf", type=str, help="name of output pdf file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # load dataframe
    df = pd.read_csv(args.df)

    # select secondary and primary reads
    sdf = df[df.sec]
    pdf = df[(~df.sec) & (~df.suppl)].set_index("read", verify_integrity=True)

    # find stetches to link
    links = []
    for i, r in sdf.iterrows():

        # reference start and end on secondary read
        rs, re = r.rs, r.re

        # reference start and end on corresponding primary read
        # the offset is given by the secondary read
        p = pdf.loc[r.read]
        ps = p.rs - p.qs + r.qs
        pe = ps + r.ref_len

        # export coordinates to link
        agree = r.fwd == p.fwd
        x = (ps, pe)
        y = (rs, re) if agree else (re, rs)
        links.append((x, y, agree))

    # Perform the plot
    fig, axs = plt.subplots(
        2, 1, figsize=(9, 11), gridspec_kw={"height_ratios": [1, 0.1]}, sharex=True
    )
    ax = axs[0]
    for x, y, c in links:
        ax.plot(x, y, alpha=0.1, color="k" if c else "r", rasterized=True)

    # legend
    legend_elements = [
        Line2D([0], [0], color=c, label=k, alpha=0.3)
        for k, c in zip(["fwd", "rev"], ["k", "r"])
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    ax.axis("equal")
    ax.set_xlabel("primary read location (bp)")
    ax.set_ylabel("secondary read location (bp)")

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
    ax.set_ylabel("n. secondary reads")

    plt.tight_layout()
    plt.savefig(args.pdf)
    plt.close(fig)
