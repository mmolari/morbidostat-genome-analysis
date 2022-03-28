# %%
import numpy as np
import pathlib as pth
import matplotlib.pyplot as plt
import re

from matplotlib.ticker import MultipleLocator


try:
    from pileupplots_utils import *
except:
    from .pileupplots_utils import *


# %%

if __name__ == "__main__":

    # argument parsing
    parser = argparser()
    args = parser.parse_args()
    data_path = pth.Path(args.vial_fld)
    fig_path = pth.Path(args.fig_fld)
    fig_path.mkdir(exist_ok=True)

    # override show and save functions
    show = show_function(args.show)
    savefig = savefig_function(fig_path)

    # %%

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing insertion plots for vial {vial}")

    # load pileups and reference genome
    insertions = load_insertions(data_path)
    ref_genome = load_reference_genome(data_path)

    # length of insertions in a dictionary { length -> n. insertions}
    ins_Ls = {k: insertion_len(ins) for k, ins in insertions.items()}

    # density of number of insertions, evaluated every `step` basepairs.
    # both as number of insertions and length of sequence inserted.
    # for each timepoints, objects are matrices whose rows correspond to
    # intervals of length `step` in the genome, and the two columns are
    # for forward and reverse reads.
    step = 5000
    L_tot = len(ref_genome)
    n_ins_dens, len_ins_dens = {}, {}
    for t, ins in insertions.items():
        n_ins_dens[t], len_ins_dens[t] = insertion_density(ins, step, L_tot)

    # create dictionary of colors for plotting
    colors = color_dict(insertions)

    # pick a threshold length above which only "n_threshold" insertions are
    # present in the last timepoint.
    n_selected = 25
    ins_len_threshold = pick_insertion_threshold_length(ins_Ls, n_selected)

    # %%
    # plot the distribution of insertion lengths and highligh the position of the longest `n_selected` insertions

    fig, axs = plt.subplot_mosaic("""AABBBBBB""", figsize=(20, 4))
    ax = axs["A"]
    max_L = max([max(Ls.keys()) for t, Ls in ins_Ls.items()])
    bins = np.logspace(0, np.log10(max_L) + 0.05, 40)
    for t, Ls in ins_Ls.items():
        l = list(Ls.keys())
        n = list(Ls.values())
        ax.hist(
            l,
            weights=n,
            label=t,
            bins=bins,
            histtype="step",
            color=colors[t],
            cumulative=-1,
        )
    ax.axvline(ins_len_threshold, ls=":", color="gray")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend(title="time")
    ax.set_xlabel("insertion length threshold [bp]")
    ax.set_ylabel("n. insertions with length < threshold")
    ax.set_title("insertion length distribution")

    ax = axs["B"]
    yticks, ylab = [], []
    for nt, t in enumerate(insertions):
        Is = insertions[t]
        x = []
        for pos, I in Is.items():
            for seq, nums in I.items():
                if len(seq) < ins_len_threshold:
                    continue
                x += [pos] * nums.sum()

        y = -np.ones_like(x) * nt
        yticks.append(-nt)
        ylab.append(t)
        ax.scatter(x, y, alpha=0.5, color=colors[t])
    ax.xaxis.set_minor_locator(MultipleLocator(100000))
    ax.grid(axis="x", which="both", alpha=0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylab)
    ax.set_title(f"location of insertions with length >= {ins_len_threshold} bp")
    ax.set_ylabel("timepoints")
    ax.set_xlabel("genome position [bp]")

    plt.tight_layout()
    savefig("insertions_length_distr.pdf")
    show()

    # %%
    # plot the density of insertions along the genome for subsequent
    # timepoints (n. insertions for forward and reverse)

    # number of timepoints
    NT = len(insertions)

    fig, axs = plt.subplots(NT, 1, figsize=(20, NT * 3), sharex=True, sharey=False)
    for t_idx, t in enumerate(insertions):
        ax = axs[t_idx]
        nins = n_ins_dens[t]
        x = np.arange(nins.shape[0]) * step
        ax.plot(x, nins[:, 0], label="fwd")
        ax.plot(x, nins[:, 1], label="rev")
        ax.set_title(f"time = {t}")

        ax.xaxis.set_minor_locator(MultipleLocator(100000))
        ax.grid(axis="both", which="both", alpha=0.5)
        ax.set_xlim(x[0], x[-1])
    axs[0].legend()
    axs[-1].set_xlabel("genome position [bp]")
    fig.supylabel("insertion density (n. insertions per bp)")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.12, left=0.05)
    savefig("insertion_density_1.pdf")
    show()

    # %%
    # plot the density of insertions along the genome for subsequent
    # timepoints (n. insertions and n. nucleotides inserted)

    # number of timepoints
    NT = len(insertions)

    fig, axs = plt.subplots(NT, 1, figsize=(20, NT * 3), sharex=True, sharey=False)
    # taxs = []
    for t_idx, t in enumerate(insertions):
        ax = axs[t_idx]
        nins, lins = n_ins_dens[t], len_ins_dens[t]
        x = np.arange(nins.shape[0]) * step
        (leg_a,) = ax.plot(x, nins.sum(axis=1), "C4")

        tax = plt.twinx(ax)
        (leg_b,) = tax.plot(x, lins.sum(axis=1), "C2")

        # ax.set_yscale("log")
        ax.set_title(f"time = {t}")
        # ax.set_xlabel("insertion length threshold [bp]")
        # ax.set_ylabel("n. insertions with length < threshold")
        ax.xaxis.set_minor_locator(MultipleLocator(100000))
        ax.grid(axis="x", which="both", alpha=0.5)
        ax.set_xlim(x[0], x[-1])

        # taxs.append(tax)

    # m = min([tax.get_ylim()[0] for tax in taxs])
    # M = max([tax.get_ylim()[1] for tax in taxs])
    # for tax in taxs:
    #     tax.set_ylim(m, M)

    axs[0].legend([leg_a, leg_b], ["n. insertions", "len. insertions"])
    axs[-1].set_xlabel("genome position [bp]")
    fig.supylabel("insertion density (n. insertions per bp)")
    fig.text(
        s="insertion density (length insertions per bp)",
        x=0.97,
        y=0.5,
        rotation="vertical",
        fontsize="large",
        verticalalignment="center",
    )
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.12, left=0.05, right=0.95)
    savefig("insertion_density_2.pdf")
    show()
