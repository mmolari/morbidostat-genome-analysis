import pysam
import functools
import pathlib
import copy

import numpy as np

from collections import defaultdict


# def pysam_open(func):
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         assert "sam" in kwargs, "function must be called with sam in kwargs"
#         sam = kwargs["sam"]
#         if isinstance(sam, str):
#             kw_copy = copy.deepcopy(kwargs)
#             with pysam.Samfile(sam) as sam_handle:
#                 kw_copy["sam"] = sam_handle
#                 return func(*args, **kw_copy)
#         else:
#             return func(*args, **kwargs)

#     return wrapper

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


# Note: the data structure for saving the pileup is a tensor, whose indices are:
# [fwd/rev, bp, L] -->  count
def ac_array(length):
    return np.zeros((2, 6, length), dtype=int)


# Note: the data structure for inserts is a nested dict with:
# position --> string  -->         count
#  (dict)      (dict)      (count fwd / count rev)
def insertion_datastruct():
    return defaultdict(lambda: defaultdict(lambda: np.zeros(2, int)))


# Note: the data structure for saving clip points is a dictionary:
# position -->  count
#  (dict)   (count fwd / count rev / tot fwd / tot rev)
def clip_datastructure():
    return defaultdict(lambda: np.zeros(4, int))
