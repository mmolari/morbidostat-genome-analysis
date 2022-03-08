# Example script

This document contains the list of commands used for the analysis of the test dataset.

Genomic analysis(Genomics):
* Identification,
* measurement,
* comparison of genomic features.

What can be measured: DNA sequence, structural variation, gene expression, annotation. Indel frequency.

[Description strongly based on the nature portfolio entry on genome analysis](https://www.nature.com/subjects/genomic-analysis)

### Description of workflow
1. Read mapping against reference genomes
2. BAM creation (merge sort and mark duplicates)
3. Indexing BAM file (Index with all duplicates marked)
4. Variant calling (Find SNPs and small indels)
5. CNV calling (Copy number changes) etc.

### SAM files
Stores seequence data in a series of tab delimetd columns. They are human readable. Generated from fastq reads and reference genome. Assign all sequences to a position with respect to the reference.

### BAM files


### CIGARS
Give additional information on the alignment of the reference and read sequence.

Gives position where alignment starts (POS: 5).

Describes alignment, is base present in both reference and read or only in one.

Example:

Pos.        1   2   3   4   5   6   7   8
Reference.  C   C   A       T   C   G   T
Read.           C   A   G   T       C   T

**CIGAR:**

POS: 2 (Alignment starts at pos 2)
CIGAR: 2M1I1M1D2M

This means: Two matches, then one insertion, then 1 match, then one deletion, at the end there are two matches (GC counts ase a match, since they still align).

**SAM and BAM files contain the same data, just in a different form**

### minimap2

```
minimap2 -xda map-ont vial1_assembly.fa reads.fq > alignment.sam
```



### SAMtools

### Pile up

### Richard scripts