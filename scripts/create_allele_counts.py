# Script adapted from https://github.com/neherlab/SVVC/blob/master/src/create_allele_counts.py
# original author: Richard Neher

import numpy as np

nuc_alpha = np.array(["A", "C", "G", "T", "-", "N"], dtype="S1")


def sam_to_allele_counts(
    sam_fname,
    qual_min=30,
    max_reads=-1,
    max_isize=700,
    VERBOSE=0,
):
    """
    calculates the allele counts for a set of mapped reads
    parameters:
    sam_fname   --   sam or bam file with mapped reads
    max_isize   --   maximal insert sizes to consider. this can be used to remove artifactual mappings
    qual_min    --   Ignore bases with quality less than qmin
    """
    import pysam
    from collections import defaultdict

    alpha = nuc_alpha

    def ac_array(length):
        return np.zeros((2, 6, length), dtype=int)

    # Note: the data structure for inserts is a nested dict with:
    # position --> string  -->         count
    #  (dict)      (dict)      (count fwd / count rev)
    def insertion_data_structure():
        return defaultdict(lambda: defaultdict(lambda: np.zeros(2, int)))

    # Open BAM or SAM file
    with pysam.Samfile(sam_fname) as samfile:
        ac = []
        refs = {}
        for nref in range(samfile.nreferences):
            if VERBOSE:
                print(
                    (
                        "allocating for:",
                        samfile.getrname(nref),
                        "length:",
                        samfile.lengths[nref],
                    )
                )
            refs[nref] = samfile.getrname(nref)
            ac.append(
                (
                    samfile.getrname(nref),
                    ac_array(samfile.lengths[nref]),
                    insertion_data_structure(),
                )
            )

        # Iterate over single reads
        for i, read in enumerate(samfile):
            # Max number of reads
            if i == max_reads:
                if VERBOSE >= 2:
                    print(("Max reads reached:", max_reads))
                break

            if (
                read.is_unmapped
                or np.abs(read.isize) > max_isize
                or read.is_secondary
                or read.is_supplementary
            ):
                continue

            # Print output
            if (VERBOSE > 2) and (not ((i + 1) % 10000)):
                print((i + 1))

            # Read CIGARs (they should be clean by now)
            counts = ac[read.rname][1][int(read.is_reverse)]
            insertion = ac[read.rname][2]

            seq = np.frombuffer(read.seq.encode(), "S1")
            qual = np.frombuffer(read.qual.encode(), np.int8) - 33
            not_primer = np.ones_like(seq, "bool")
            pos = read.pos

            # if pos+len(seq)>7267:
            # 	import ipdb;ipdb.set_trace()
            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):
                if block_type == 4:  # softclip
                    seq = seq[block_len:]
                    qual = qual[block_len:]
                    # not the difference here: the reported position starts after the softclip. hence the not_primer is already correct
                    not_primer = not_primer[:-block_len]
                    continue
                if block_type == 5:  # hard clip
                    continue

                # Check for pos: it should never exceed the length of the fragment
                #                if (block_type in [0, 1, 2]) and (pos >= length):
                #                    raise ValueError('Pos exceeded the length of the fragment')

                # Inline block
                if block_type == 0:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    not_primerb = not_primer[:block_len]
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = (
                            (seqb == a) & (qualb >= qual_min) & (not_primerb)
                        ).nonzero()[0]
                        if len(posa):
                            counts[j, pos + posa] += 1

                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        not_primer = not_primer[block_len:]
                        pos += block_len

                # Deletion
                elif block_type == 2:
                    # Increment gap counts
                    counts[4, pos : pos + block_len] += 1
                    # Chop off pos, but not sequence
                    pos += block_len

                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    not_primerb = not_primer[:block_len]
                    # Accept only high-quality inserts
                    if (qualb >= qual_min).all():
                        insertion[pos][seqb.tobytes().decode()][
                            int(read.is_reverse)
                        ] += 1

                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        not_primer = not_primer[block_len:]

                # Other types of cigar?
                else:
                    if VERBOSE > 2:
                        print(("unrecognized CIGAR type:", read.cigarstring))
                    # raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return ac


def dump_allele_counts(dirname, ac, suffix=""):
    import pickle, gzip, os

    dirname = dirname.rstrip("/") + "/"
    if not os.path.isdir(dirname):
        print(("creating directory", dirname))
        try:
            os.mkdir(dirname)
        except:
            raise "creating directory failed"

    for refname, ac_array, insertions in ac:
        print(refname)
        np.savez_compressed(dirname + "allele_counts" + suffix + ".npz", ac_array)
        with gzip.open(dirname + "insertions" + suffix + ".pkl.gz", "w") as outfile:
            pickle.dump({k: dict(v) for k, v in insertions.items()}, outfile)


def load_allele_counts(dirname, suffix="", allCounts=False):
    import pickle, gzip, glob

    dirname = dirname.rstrip("/") + "/"
    tmp_ac = {}
    if allCounts:
        ac_flist = glob.glob(dirname + "*allele_counts" + suffix + ".npz")
    else:
        ac_flist = glob.glob(dirname + "allele_counts" + suffix + ".npz")
    for fname in ac_flist:
        # print("reading",fname)
        tmp = "_allele_counts" + suffix + ".npz"
        refname = fname.split("/")[-1][: -len(tmp)]
        tmp_ac[refname] = list(np.load(fname).items())[0][1]

    ins_flist = glob.glob(dirname + "*insertions" + suffix + ".pkl.gz")
    tmp_ins = {}
    for fname in ins_flist:
        # print("reading",fname)
        tmp = "_insertions" + suffix + ".pkl.gz"
        refname = fname.split("/")[-1][: -len(tmp)]
        with gzip.open(fname) as fh:
            tmp_ins[refname] = pickle.load(fh)

    ac = []
    for refname in tmp_ac:
        ac.append((refname, tmp_ac[refname].sum(axis=0).sum(axis=0)))

    # ins = []
    # for refname in tmp_ins:
    #     ins.append((refname, tmp_ins[refname]))

    return ac, tmp_ins


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="create allele counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam_file", help="bam file to pile up")
    parser.add_argument("--out_dir", help="directory to save results")
    parser.add_argument("--qual_min", help="minimum acceptable read quality", type=int)
    parser.add_argument(
        "--max_insertionsize", help="maximum insertion size saved", type=int
    )

    args = parser.parse_args()

    ac = sam_to_allele_counts(
        args.bam_file,
        qual_min=args.qual_min,
        VERBOSE=3,
        max_isize=args.max_insertionsize,
    )

    ac_renamed = []
    for refname, counts, insertions in ac:
        ac_renamed.append((refname.replace("/", "_"), counts, insertions))
    ac = ac_renamed

    dump_allele_counts(args.out_dir, ac)
