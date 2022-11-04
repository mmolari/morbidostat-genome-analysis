# %%
import pathlib as pth

try:
    from utils.extract_stats_utils import *
except:
    from .utils.extract_stats_utils import *

# %%

if __name__ == "__main__":

    parser = argparser()
    parser.description = """
    Evaluates the frequency of reads that agree with the reference genome,
    for every timepoint. The reference genome is the one corresponding to
    the first timepoint. Expects an input folder containing data for one
    vial, organized in timepoints. The first vial must contain a reference genome.
    It produces a stat_table_reference_freq.pkl.gz file in the same folder.

    Statistic: agreement with reference genome.
    Ntrue: n. reads that have the same value as the reference genome.
    Ntot: total number of reads per site (no 'N' or '-').
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

    vprint(f"extracting reference frequency for {input_fld}")

    # load pileups and reference genome
    vprint(f"load pileup and reference genome")
    pileups = load_pileups(input_fld)
    ref_genome = load_reference_genome(input_fld)

    # evaluate number of times a nucleotide is consensus and total number of nucleotides:
    vprint(f"evaluate consensus frequency")
    ncons, ntot = {}, {}
    for k, pp in pileups.items():
        ncons[k], ntot[k] = consensus_count(pp, ref_genome)

    vprint(f"create StatsTable")
    st = StatsTable.from_counts(
        stat_name="reference_freq",
        Ntrue=ncons,
        Ntot=ntot,
    )

    vprint(f"save StatsTable")
    st.save(fld=output_fld)

    vprint(f"done")
