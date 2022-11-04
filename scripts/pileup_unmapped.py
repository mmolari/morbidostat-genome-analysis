import numpy as np
import pysam
import argparse
import json

from utils.pileup_utils import ref_info


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=str, help="sorted bam file")
    parser.add_argument("--json_out", type=str, help="output json file")
    return parser.parse_args()


# %%
if __name__ == "__main__":

    args = parse_args()

    # list of unmapped reads
    unmapped = []

    with pysam.Samfile(args.bam) as sf:

        # capture refs in reads
        refs = ref_info(sf)
        assert len(refs) == 1, f"Only one reference allowed for pileup: {refs}"
        ref = refs[0]

        for i, read in enumerate(sf):

            # capture unmapped reads
            if read.is_unmapped:
                unmapped.append(
                    {
                        "read": read.query_name,
                        "flag": read.flag,
                        "len": len(read.seq),
                        "avg. qscore": np.mean(read.query_qualities),
                        "seq": read.seq,
                        "qscore": read.qual,
                    }
                )

    with open(args.json_out, "w") as f:
        json.dump(unmapped, f)
