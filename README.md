# morbidostat-genome-analysis

Code for the analysis of morbidostat data. The analysis consists for now of the following steps:

- map reads against assembled genome from the first timepoint
- produce a pileup of the mapped reads, that can be explored with IGV

## content of the repository

- The folder `notes` contains notes and descrioptions of the steps implemented in the pipeline. To be used as guides for pipeline implementation.

## test dataset

A test dataset is available to test the steps of the pipeline. This can be downloaded from scicore cluster using:
```
username@login.scicore.unibas.ch:/
```
Substituting `username`

The test dataset has been extracted from the morbidostat run labeled `2022-02-08-RT_test`. It contains 2 vials, each one having 3 timepoints. For the first timepoint the annotated genome is included. For the rest of the timepoints, only fastq reads are provided, and these files contain only a subset of all the reads (first 10'000).