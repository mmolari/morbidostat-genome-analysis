import pysam
import functools
import pathlib
import copy

import numpy as np

from collections import defaultdict


nuc_alpha = np.array(["A", "C", "G", "T", "-", "N"], dtype="S1")


def ref_info(samfile):
    Nrefs = samfile.nreferences
    refs = [
        {"name": samfile.getrname(n), "len": samfile.lengths[n]} for n in range(Nrefs)
    ]
    return refs


def read_info_dict(read):
    cigar = read.cigartuples
    Hstart = 0
    if cigar[0][0] == 5:
        Hstart = cigar[0][1]
    res = {
        "read": read.query_name,
        "flag": read.flag,
        "fwd": read.is_forward,
        "ref_len": read.reference_length,
        "qry_len": read.query_length,
        "sec": read.is_secondary,
        "suppl": read.is_supplementary,
        "rs": read.reference_start,
        "re": read.reference_end,
        "qs": read.query_alignment_start + Hstart,
        "qe": read.query_alignment_end + Hstart,
    }
    return res