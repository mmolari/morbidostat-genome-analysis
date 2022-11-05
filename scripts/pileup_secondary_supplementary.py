import pysam
import argparse

import numpy as np
import pandas as pd

from utils.pileup_utils import ref_info, read_info_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that given a sam/bam file, creates csv dataframe on info on secondary
    and supplementary mappings, with their corresponding primary one.
    """
    )
    parser.add_argument("--bam", type=str, help="sorted bam file")
    parser.add_argument("--df_out", type=str, help="output csv dataframe")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # list of unmapped reads
    non_primary_ids = []
    df = []

    with pysam.Samfile(args.bam) as sf:

        # capture refs in reads
        refs = ref_info(sf)
        assert len(refs) == 1, f"Only one reference allowed for pileup: {refs}"
        ref = refs[0]

        for read in sf.fetch():
            if read.is_secondary or read.is_supplementary:
                non_primary_ids.append(read.query_name)
        non_primary_ids = np.unique(non_primary_ids)

        for read in sf.fetch():
            if read.query_name in non_primary_ids:
                df.append(read_info_dict(read))

    if len(non_primary_ids) > 0:
        df = pd.DataFrame(df).sort_values(["read", "flag"]).reset_index(drop=True)
        df.to_csv(args.df_out, index=False)
