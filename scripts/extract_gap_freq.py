# %%
from lib2to3.pytree import NegatedPattern
import pathlib as pth

try:
    from extract_stats_utils import *
except:
    from .extract_stats_utils import *

# %%

if __name__ == "__main__":

    parser = argparser()
    parser.description = """
    Evaluates the frequency of gaps per position in the reads for every timepoint.
    Expects an input folder containing data for one vial, organized in timepoints.
    It produces a stat_table_gap_freq.pkl.gz file in the same folder.

    Statistic: gap frequency.
    Ntrue: n. reads that have a gap in a given position.
    Ntot: total number of reads on the position (including '-' and 'N').
    """
    args = parser.parse_args()
    input_fld = pth.Path(args.vial_fld)
    output_fld = pth.Path(args.output_fld)
    output_fld.mkdir(exist_ok=True)

    if args.verbose:
        vprint = print
    else:
        vprint = lambda x: None

    # %%
    # input_fld = pth.Path("../results/2022-02-08_RT_test/vial_04/")
    # vprint = print

    vprint(f"extracting gap frequency for {input_fld}")

    # load pileups and reference genome
    vprint(f"load pileup")
    pileups = load_pileups(input_fld)

    # evaluate number of times a nucleotide is consensus and total number of nucleotides:
    vprint(f"evaluate gap frequency")

    # get gap counts
    Ngap, Ntot = {}, {}
    for k, pp in pileups.items():
        Ngap[k], Ntot[k] = gap_count(pp)

    vprint(f"create StatsTable")
    st = StatsTable.from_counts(
        stat_name="gap_freq",
        Ntrue=Ngap,
        Ntot=Ntot,
    )

    vprint(f"save StatsTable")
    st.save(fld=output_fld)

    vprint(f"done")
