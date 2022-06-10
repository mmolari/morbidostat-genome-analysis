# morbidostat-genome-analysis

Code for the analysis of morbidostat data. The analysis consists for now of the following steps:

- map reads against assembled genome from the first timepoint
- produce a pileup of the mapped reads, that can be explored with IGV

## content of the repository

- The folder `notes` contains notes and descrioptions of the steps implemented in the pipeline. To be used as guides for pipeline implementation.
- The `scripts` folder contains scripts used by the various workflows.

## test dataset

A test dataset is available to test the steps of the pipeline. This can be downloaded from scicore cluster using:
```
scp -r username@login.scicore.unibas.ch:/scicore/home/neher/GROUP/data/2022_nanopore_sequencing/experiments/test_dataset .
```
Substituting `username` with the appropriate username.

The test dataset has been extracted from the morbidostat run labeled `2022-02-08-RT_test`. It contains 2 vials (vial 2 and 3), each one having 3 timepoints (timepoint 1, 2 and 5). For timepoint 1 the annotated genome is included. For the rest of the timepoints, only fastq reads are provided in the `reads.fastq.gz` files. These files contain only a subset of all the reads (first 10'000 lines of the original fastq file). This corresponds to ~3% of all the reads, and these files are roughly 20-30 Mb each.

## workflow: building a pileup of mapped reads

The workflow `pileup.nf` builds a pileup of the reads. For every vial, it maps reads for any timepoints to a reference genome, which is the one assembled from the initial timepoint for each vial. It generates the sorted bam file of mappings, together with a pileup matrix and a list of insertions. Usage is as follows:

```bash
nextflow run pileup.nf \
    -profile cluster \
    --input_fld test_dataset \
    --time_beg 1 \
    --qual_min 15 \
    --max_insertion_size 1000  \
    -resume
```

The input parameters are:
- `profile` : if `cluster` is specified, then SLURM execution is activated. Otherwise local execution is used.
- `input_fld` : folder containing the input reads. It must have the nested `vial_n/time_n` structure of archived morbidostat experiment data.
- `time_beg` : the label corresponding to the first timepoint, whose assembled genome is used as reference (e.g. `1` for timepoint_1`)
- `qual_min` and `max_insertion_size` : parameters of the script used to build the pileup. Only reads with quality higher than the threshold are used, and only insertions shorter than the threshold are considered.

results are saved in the `results` folder, with a structure that mirrors the `vial_n/timepoint_n` structure of the input folder. Each of these subfolders will contain the following files:
- `reads.sorted.bam` and `reads.sorted.bam.bai` : sorted `bam` file (and corresponding index) containing the mapping of the reads against the reference genome.
- `pileup/allele_counts.npz` : pileup of the reads. This is a numpy tensor with dimension (2,6,L) corresponding to (1) forward-reverse reads, (2) allele `["A", "C", "G", "T", "-", "N"]` and (3) position.
- `pileup/insertions.pkl.gz` : a nested dictionary of insertions, saved in pickle format and compressed with gzip. The structure is `position -> sequence -> [n. forward reads, n. reverse reads]`.
- `ref_genome.fa` and `ref_genome.gbk` : symlink to the reference genome used for mapping the reads (same vial, first timepoint), both in genbank and fasta format. 


## workflow: extract statistics

The workflow `extract_stats.nf` extracts relevant statistics from the pileups and insertion lists. These include:

- frequency of reads being equal to the reference genome
- frequency of gaps

These are saved as `stats/stats_table_{statname}.pkl.gz` files containing pandas dataframes.

The workflow can be run by executing:

```
nextflow run extract_stats.nf \
    -profile cluster \
    --input_fld results/test_dataset \
    -resume
```

where `--input_fld` is the subfolder of `results` corresponding to the dataset of interest, containing the pileups for different vials and different timepoints.

## workflow: summary plots with pileup analysis

The workflow `plots.nf` executes scripts to analyze the results of the pileup. The workflow can be launched with:

```bash
nextflow run plots.nf \
    -profile cluster \
    --input_fld results/test_dataset \
    -resume
```

The `--input_fld` flag is used to specify the folder in which the results of the previous workflow are saved.
The output figures are saved in `figures/subfolder` where the subfolder has the name of the run, and they are separated by vial.
The executed scripts are:
- `pileupplots_coverage.py`
- `pileupplots_consensus_frequency.py`
- `pileupplots_gaps.py`
- `pileupplots_insertions.py`

For each of these scripts, usage is:

```bash
python3 scripts/pileupplots_consensus_frequency.py \
    --vial_fld results/2022-02-08_RT_test/vial_02 \
    --fig_fld     figs/2022-02-08_RT_test/vial_02 \
    --show
```

Where:
- `--vial_fld` is the folder containing the results for a single vial (pileup, reference genome...)
- `--fig_fld` is the output folder in which to save figures (must exist)
- `--show` if specified figures are visualized with `plt.show()`, and otherwise they are simply saved.