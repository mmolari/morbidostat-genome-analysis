import pysam
import argparse
import sys
import gzip

import numpy as np
import pandas as pd

from Bio import SeqIO, SeqRecord, Seq

from utils.pileup_utils import ref_info


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that given a sam/bam file, creates a .fastq.gz file with unmapped reads
    and a .csv dataframe with stats on each read (name, length, avg_quality, flag).
    """
    )
    parser.add_argument("--bam", type=str, help="sorted bam file")
    parser.add_argument("--df_out", type=str, help="output csv dataframe")
    parser.add_argument("--fastq_out", type=str, help="output fastq.gz file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # list of unmapped reads
    unmapped_df = []
    unmapped_records = []

    with pysam.Samfile(args.bam) as sf:

        # capture refs in reads
        refs = ref_info(sf)
        assert len(refs) == 1, f"Only one reference allowed for pileup: {refs}"
        ref = refs[0]

        for i, read in enumerate(sf):

            # capture only unmapped reads
            if not read.is_unmapped:
                continue

            unmapped_df.append(
                {
                    "read": read.query_name,
                    "len": len(read.seq),
                    "avg. qscore": np.mean(read.query_qualities),
                    "flag": read.flag,
                }
            )

            # create and append fastq-compatible record
            seq = Seq.Seq(read.seq)
            rec = SeqRecord.SeqRecord(seq, id=read.query_name, description="", name="")
            rec.letter_annotations["phred_quality"] = read.query_qualities
            unmapped_records.append(rec)

    if len(unmapped_df) == 0:
        sys.exit(0)

    # create and export dataframe
    df = pd.DataFrame(unmapped_df).sort_values("len", ascending=False)
    df.to_csv(args.df_out, index=False)

    # create and export fastq.gz file
    with gzip.open(args.fastq_out, "wt") as f:
        SeqIO.write(unmapped_records, f, "fastq")
