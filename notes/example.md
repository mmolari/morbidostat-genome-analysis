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
Binary files (= Non human readable). Contain the same information as a SAM file but have a smaller size and are therefore in many cases preferrable to SAM files. Can have sorted and non sorted reads.

### CIGARS
Give additional information on the alignment of the reference and read sequence.

Gives position where alignment starts (POS: 5).

Describes alignment, is base present in both reference and read or only in one.

Example:
```
Pos.        1   2   3   4   5   6   7   8
Reference.  C   C   A       T   C   G   T
Read.           C   A   G   T       C   T
```
**CIGAR:**

POS: 2 (Alignment starts at pos 2)
CIGAR: 2M1I1M1D2M

This means: Two matches, then one insertion, then 1 match, then one deletion, at the end there are two matches (GC counts ase a match, since they still align).

**SAM and BAM files contain the same data, just in a different form**

### minimap2

Maps reads to the reference genome (for us usually the one from the first timepoint in the experiment).

```
minimap2 -a -x -t map-ont vial1_assembly.fa reads.fq > alignment.sam


-a for alignment
-d for indexing
-t number of threads
-x for nanopore reads, same as map-ont
```
`May want to add -L tag for ultra long reads`

minmap outputs up to 5 other alignments (pseudogenes). The output is in the SAM format.

### SAMtools

1. Transforms the mapped reads into a sorted BAM file.

```
samtools sort -@ threads reads.sam > reads.sorted.bam
```

2. Creating an index for the sorted BAM file.

```
samtools index -@ threads reads.sorted.bam
```

3. (Option: Create pileup with samtools)

### Pile up (Richards scripts)


