# Script adapted from https://github.com/neherlab/SVVC/blob/master/src/create_allele_counts.py
# original author: Richard Neher

# %%
import numpy as np

nuc_alpha = np.array(["A", "C", "G", "T", "-", "N"], dtype="S1")


def sam_to_allele_counts(
    sam_fname,
    qual_min=30,
    clip_minL=100,
    VERBOSE=0,
):
    """
    calculates the allele counts for a set of mapped reads
    parameters:
    sam_fname   --  sam or bam file with mapped reads
    qual_min    --  Ignore bases with quality less than qmin
                    and insertions with less than 70% reads meeting this threshold
    clip_minL   --  ignore clipped reads with less than this threshold quality
    """
    import pysam
    from collections import defaultdict

    alpha = nuc_alpha

    def ac_array(length):
        return np.zeros((2, 6, length), dtype=int)

    # Note: the data structure for inserts is a nested dict with:
    # position --> string  -->         count
    #  (dict)      (dict)      (count fwd / count rev)
    def insertion_datastruct():
        return defaultdict(lambda: defaultdict(lambda: np.zeros(2, int)))

    # Note: the data structure for saving clip points is a dictionary:
    # position -->  count
    #  (dict)   (count fwd / count rev)
    def clip_datastructure():
        return defaultdict(lambda: np.zeros(2, int))

    # Open BAM or SAM file
    with pysam.Samfile(sam_fname) as samfile:
        ac = []

        # allocate space
        for nref in range(samfile.nreferences):
            name = samfile.getrname(nref)
            L = samfile.lengths[nref]
            if VERBOSE:
                print(("allocating for:", name, "length:", L))
            ac.append((name, ac_array(L), insertion_datastruct(), clip_datastructure()))

        # Iterate over single reads
        for i, read in enumerate(samfile):

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Read CIGARs (they should be clean by now)
            counts = ac[read.rname][1][int(read.is_reverse)]
            insertion = ac[read.rname][2]
            clips = ac[read.rname][3]

            seq = np.frombuffer(read.seq.encode(), "S1")
            qual = np.frombuffer(read.qual.encode(), np.int8) - 33
            pos = read.pos

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):

                # Print output
                if (VERBOSE > 2) and (not ((ic + 1) % 10000)):
                    print(f"read n. {ic + 1}")

                if block_type == 4:  # softclip
                    if block_len >= clip_minL:
                        clips[pos][int(read.is_reverse)] += 1
                    seq = seq[block_len:]
                    qual = qual[block_len:]
                    continue
                if block_type == 5:  # hard clip
                    if block_len >= clip_minL:
                        clips[pos][int(read.is_reverse)] += 1
                    continue
                if block_type == 0:  # Inline block
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                        if len(posa):
                            counts[j, pos + posa] += 1

                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        pos += block_len

                elif block_type == 2:  # Deletion
                    # Increment gap counts
                    counts[4, pos : pos + block_len] += 1
                    # Chop off pos, but not sequence
                    pos += block_len

                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:  # insertion
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    # Accept only high-quality inserts
                    if (qualb >= qual_min).mean() > 0.7:
                        insertion[pos][seqb.tobytes().decode()][
                            int(read.is_reverse)
                        ] += 1

                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]

                # Other types of cigar?
                else:
                    if VERBOSE > 2:
                        print(("unrecognized CIGAR type:", read.cigarstring))

    return ac


# %%
def dump_allele_counts(dirname, ac, suffix=""):
    import pickle, gzip, os

    dirname = dirname.rstrip("/") + "/"
    if not os.path.isdir(dirname):
        print(("creating directory", dirname))
        try:
            os.mkdir(dirname)
        except:
            raise "creating directory failed"

    for refname, ac_array, insertions, clips in ac:
        print(refname)
        np.savez_compressed(dirname + "allele_counts" + suffix + ".npz", ac_array)
        with gzip.open(dirname + "insertions" + suffix + ".pkl.gz", "w") as outfile:
            pickle.dump({k: dict(v) for k, v in insertions.items()}, outfile)
        with gzip.open(dirname + "clips" + suffix + ".pkl.gz", "w") as outfile:
            pickle.dump(dict(clips), outfile)


# %%
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="create allele counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam_file", help="bam file to pile up")
    parser.add_argument("--out_dir", help="directory to save results")
    parser.add_argument("--qual_min", help="minimum acceptable read quality", type=int)
    parser.add_argument("--clip_minL", help="minimum acceptable read quality", type=int)

    args = parser.parse_args()

    ac = sam_to_allele_counts(
        args.bam_file, qual_min=args.qual_min, clip_minL=args.clip_minL, VERBOSE=3
    )

    ac_renamed = []
    for refname, counts, insertions, clips in ac:
        ac_renamed.append((refname.replace("/", "_"), counts, insertions, clips))
    ac = ac_renamed

    dump_allele_counts(args.out_dir, ac)
