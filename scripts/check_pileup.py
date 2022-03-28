# %%

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.stats import entropy


# %%

# vial number and timepoint
vn, tn = 2, 5
subfld = pathlib.Path(f"../results/2022-02-08_RT_test/vial_{vn:02d}/time_{tn:d}")

# load pileup file
pileup = np.load(subfld / "pileup/allele_counts.npz")["arr_0"]

# load fasta file with genome
fasta_file = subfld / "ref_genome.fa"
with open(fasta_file, "r") as f:
    ref_genome = SeqIO.read(f, format="fasta")

# %%
def coverage_fwd_rev(pp):
    """Extract coverage (only ACGT, no gaps and N)"""
    fwd, rev = pp[0], pp[1]
    fwd_cov = np.sum(fwd[:4, :], axis=0)
    rev_cov = np.sum(rev[:4, :], axis=0)
    return fwd_cov, rev_cov


fwd_cov, rev_cov = coverage_fwd_rev(pileup)


def sum_fwd_rev(pp):
    """Sum forward and reverse reads, remove one dimension"""
    return np.sum(pp, axis=0)


tot_pp = sum_fwd_rev(pileup)


def pileup_entropy(tot_pileup):
    """Evaluates the entropy of the pileup (dimension 4xL)"""
    return entropy(tot_pileup[:6, :], axis=0, base=6)


entr = pileup_entropy(tot_pp)


nuc_alpha = ["A", "C", "G", "T", "-", "N"]
nuc_idx = {n: i for i, n in enumerate(nuc_alpha)}


def seq_to_idx(seq):
    nseq = np.array(seq)
    return np.array([nuc_idx[x] for x in nseq])


refgen_idx = seq_to_idx(ref_genome)


def fraction_consensus(pileup, refgen_idx):
    tot_pileup = np.sum(pileup, axis=0)
    tot_pileup = tot_pileup[:4, :]
    tot = np.sum(tot_pileup, axis=0)
    cons = np.choose(refgen_idx, tot_pileup)
    return cons / tot


f_cons = fraction_consensus(pileup, refgen_idx)

# %%
fig, axs = plt.subplot_mosaic("ABBB", figsize=(15, 4))


ax = axs["A"]
bins = np.arange(300)
kw = {"bins": bins, "histtype": "step"}
ax.hist(fwd_cov, label="fwd", **kw)
ax.hist(rev_cov, label="rev", **kw)
ax.hist(fwd_cov + rev_cov, label="tot", **kw)
ax.set_title(f"vial {vn}, timepoint {tn}")
ax.legend()
ax.set_xlabel("coverage per site")
ax.set_ylabel("n. sites")

ax = axs["B"]
L = len(fwd_cov)
x = np.arange(L)
bins = x[::1000]
kw = {"bins": bins, "histtype": "step"}
ax.hist(x, weights=fwd_cov / 1000, **kw)
ax.hist(x, weights=rev_cov / 1000, **kw)
ax.hist(x, weights=(fwd_cov + rev_cov) / 1000, **kw)
ax.set_ylabel("coverage per site (running average over 1kbp)")
ax.set_xlabel("position on genome (bp)")

plt.tight_layout()
plt.show()

# %%

thr = 0.4

fig, axs = plt.subplot_mosaic("ABBB", figsize=(15, 4))

ax = axs["A"]
ax.hist(entr, bins=100, histtype="step")
ax.set_yscale("log")
ax.axvline(thr, ls=":", c="k")
ax.set_xlabel("entropy")
ax.set_ylabel("n. sites")
ax.set_title(f"vial {vn}, timepoint {tn}")


ax = axs["B"]
L = len(fwd_cov)
x = np.arange(L)
mask = entr >= thr
ax.hist(x[mask], bins=x[::1000], histtype="step")
ax.set_xlabel("position on genome (bp)")
ax.set_ylabel(f"n. sites with entr>{thr} per kbp")

plt.tight_layout()
plt.show()
# %%


thr = 0.4

fig, axs = plt.subplot_mosaic("ABBB", figsize=(15, 4))

ax = axs["A"]
ax.hist(f_cons, bins=100, histtype="step")
ax.set_yscale("log")
ax.axvline(thr, ls=":", c="k")
ax.set_xlabel("fraction of consensus reads")
ax.set_ylabel("n. sites")
ax.set_title(f"vial {vn}, timepoint {tn}")


ax = axs["B"]
L = len(fwd_cov)
x = np.arange(L)
mask = f_cons <= thr
# ax.hist(x[mask], bins=x[::1000], histtype="step")
ax.plot(x[mask], f_cons[mask], ".")
ax.set_xlabel("position on genome (bp)")
ax.set_ylabel(f"n. sites with fract consensus<{thr} per kbp")
# ax.set_xlim(int(1e6), int(1e6) + 1000)

plt.tight_layout()
plt.show()
# %%

# %%
