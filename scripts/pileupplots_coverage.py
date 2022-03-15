# %%
import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re


try:
    from pileupplots_utils import *
except:
    from .pileupplots_utils import *


def average_every_step(arr, step):
    """Return an average of the array every `step` positions. Also
    returns the position corresponding to each average."""
    L = arr.size
    floor = (L // step) * step
    x = (np.arange(L // step) * step) + step // 2
    subarr = arr[:floor].reshape(-1, step).mean(axis=1)
    return x, subarr


# %%


if __name__ == "__main__":

    parser = argparser()
    args = parser.parse_args()
    data_path = pth.Path(args.vial_fld)
    fig_path = pth.Path(args.fig_fld)
    fig_path.mkdir(exist_ok=True)

    # override show and save functions
    show = show_function(args.show)
    savefig = savefig_function(fig_path)

    # get vial number
    vial = re.search("vial_(\d+)/?$", str(data_path)).groups()[0]
    print(f"preparing coverage plots for vial {vial}")

    # load pileups and reference genome
    pileups = load_pileups(data_path)

    # evaluate coverages
    coverages = {k: coverage(pp) for k, pp in pileups.items()}

    # assign colors to timepoints
    colors = color_dict(pileups)
    # %%

    # coverage plot
    fig, axs = plt.subplot_mosaic("ABBBB", figsize=(20, 4))

    # histogram of coverages
    ax = axs["A"]
    maxcov = max([max(cov) for cov in coverages.values()])
    bins = np.arange(maxcov)
    cumulative_histograms(
        coverages, ax, colors, plotmeans=True, bins=bins, cumulative=True
    )
    ax.legend(loc="upper left", title="time")
    ax.set_xlabel("coverage per site")
    ax.set_ylabel("n. sites (cumulative)")

    # coverage along the genome (mean of every kbp)
    ax = axs["B"]
    ax.set_title(f"vial {vial}")
    step = 1000
    for tp, cov in coverages.items():
        cov_m = cov / cov.mean()
        x, subcov = average_every_step(cov_m, step)
        ax.plot(x, subcov, label=tp, color=colors[tp], alpha=0.8)

    # guidelines (mean and 2 std)
    cov_mat = np.vstack(list(coverages.values()))
    cov_rescaled = (cov_mat.T / cov_mat.mean(axis=1)).T
    cov_mean = cov_rescaled.mean(axis=0)
    cov_std = cov_rescaled.std(axis=0)
    step = 50000
    xm, cov_mean = average_every_step(cov_mean, step)
    xs, cov_std = average_every_step(cov_std, step)
    ax.plot(xm, cov_mean, "--", color="gray")
    ax.plot(xm, cov_mean - 2 * cov_std, ":", color="gray")
    ax.plot(xm, cov_mean + 2 * cov_std, ":", color="gray")

    ax.set_xlabel("position on the genome (bp)")
    ax.set_ylabel(f"coverage ({step} bp mean) / mean coverage")
    plt.tight_layout()
    savefig("coverage.pdf")
    show()


# %%
