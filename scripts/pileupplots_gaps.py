# %%
import numpy as np
import pathlib as pth
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

from collections import defaultdict

try:
    from pileupplots_utils import *
except:
    from .pileupplots_utils import *

# %%

# path for a single vial
# input_data_path = "../results/2022-02-08_RT_test/vial_02/"

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
    print(f"preparing frequency plots for vial {vial}")

    # load pileups and reference genome
    pileups = load_pileups(data_path)
    ref_genome = load_reference_genome(data_path)

    # evaluate consensus frequencies
    ref_genome_idxs = ref_genome_to_nucl_id(ref_genome)
    cons_freqs = {
        k: consensus_frequency(pp, ref_genome_idxs) for k, pp in pileups.items()
    }
