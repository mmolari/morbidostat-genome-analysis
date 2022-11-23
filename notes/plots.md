# Plot description

This file contains a succint descriptions of the plot and files produced by the `plot.nf` workflow. For each plotting script we give a short description of the figures and files produced. Note that the analysis is mainly based on the pileup and that it only contains reads of sufficient quality.

## Coverage

The script `plot_coverage.py` analyzes the coverage, defined as number of reads per site (excluding `-` or `N`).

- `coverage_distribution.pdf`:
    - panel 1: cumulative distribution of coverage per site for all the different timepoints. Vertical lines represent medians.
    - panel 2: coverage along the genome for the different timepoints. To compare the distribution values for the coverage are rescaled using the average coverage (so that on average the rescaled coverage is 1) and to reduce noise the average over 1kbp segments is reported. The gray dashed line represent the mean $\pm$ 2 standard deviation over a 50kbp window over all the timepoints. 
- `coverage_fwd_rev.pdf`:
    - boxplot with the distribution of forward and reverse coverage for the different timepoints.


## Consensus frequency

The script `plot_consensus_freq.py` analyzes for each site the *consensus frequency* (or better, the *reference frequency*), i.e. the fraction of times that reads mapping on a different position are the same as in the reference genome.

This script select the 99 positions that have the highest absolute consensus frequency change over time. These are the reads that are more likely to have undergone a mutation.
The absolute change is defined as the (max - min) consensus frequency over time. Timepoints that do not meet the following criteria are excluded:
- have at least 10 fwd and 10 rev reads.
- positions in which the fisher exact test for n. of consensus vs non-consensus reads restitutes a p-value > 5%. This is done to exclude reads with a strong fwd/rev disagreement.

The produced figures are:

- `consensus_freq_distribution.pdf`:
    - panel 1: distribution of consensus frequencies for the different timepoints
    - panel 2: distribution of $\Delta f = f(t_{end}) - f(t_{beg})$. Reads with minimum $\Delta f$ are likely to have undergone a mutation.
- `consensus_pvalue_distribution.pdf`:
    - joint distribution of the forward and reverse p-values, complete with marginal distribuitons. The red line represents the threshold $p_f \times p_r < 5\%$ below which positions are excluded.
- `consensus_pvalue_vs_freqdiff.pdf`: 
    - joint distribution of rank $\Delta f$ and total p-value $p_f \times p_r$, subdivided in two categories: positions that have more or less than the threshold value of 10 reads.
- `consensus_freq_trajs.pdf`:
    - each plot represents the frequency trajectory of a selected position, ranked by $\Delta f$. The title reports the position (1-based numbering) and the value of $\Delta f$. The black like represents the total frequency, and blue and orange markers report respectively forward and reverse frequencies. If a position has less than 10 read then the marker is displayed with a higher transparency.

The csv file `consensus_selected_positions.csv` has three lines per relevant positions. Positions are ordered by $\Delta f$. The three lines correspond to total, forward and reverse reads. For each line and each timepoint it contains the frequency of consensus reads `f_t` and the total number of reads `N_t`. It also contains the p-value and ranking parameter $\Delta f$.

## Gap frequency

The script `plot_gaps.py` performs a similar analysis to the one for the *consensus frequency*, but instead of the fraction of consensus reads this will consider the fraction of reads that contain a gap.

Two plots are different from consensus plots:
- `gap_vs_genome.pdf`: cumulative gap frequency over a window of 100bps, stratified by time and fwd/rev reads. Gaps that form over time should be visible as an increase in this frequency.



## Insertions

The script `plot_insertions.py` will look at insertion in the reads.

The produced figures are:
- `insertions_distributions.pdf`: the distribution of number of insertions, insertion frequency and average insertion length for the different timepoints
- `insertions_joint_distr.pdf`: for the last available timepoint plots the joint distributions of values for the forward and reverse reads for three different quantitites: number of insertions, insertion frequency and insertion length.
- `insertions_vs_genome_{L,N}.pdf`:
    - n. of insertions on a window of 100 kbp size, stratified by timepoint and fwd/rev reads.
    - same as previous, but instead of n. of reads it produces an histogram of average insertion length (L * insertion freq).
- `insertions_trajectories.pdf`: plots the insertion frequency trajectory for at different timepoints for few selected positions. The meaning of markers and colors is similar to the one for the consens frequency trajectories. Positions are selected based on the (max - min) insertion frequency over time. Only points with not less than 8 insertions are considered, and points where a fisher exact test on the number of fwd vs rev insertions returns a p-value greater than 5%.

The relevant positions are recorded in `insertions_selected_positions.csv`. This is similar to the csv file for the consensus frequency, with the addition of the `L` columns storing the average insertion length.


## Clipped reads

The script `plot_clips.py` is used to analyze soft and hard clips. In this analysis we select positions with a high number of clips.

The produced figures are:
- `clip_distributions.pdf`: distributions of n. of clips per site, clip frequency and average length of soft-clips for the different timepoints.
- `clip_joint_dist_final_time.pdf`: joint distributions of the three above quantities for forward and reverse reads, evaluated at the final timepoint.
- `clips_vs_genome.pdf`: distributions of n. of clips in forward and reverse reads over intervals spanning 5 kbps. One distribution for timepoint.
- `clips_trajectories.pdf`: trajectories of number of clips for the selected positions. Notice that this is the total number and is not normalized by the total number of reads.

The statistics for the number of clips for the relevant positions are saved in the `clips_selected_positions.csv` file, the structure is analogous to the csv file for insertions.

## Secondary/Supplementary alignments

These are all alignments for which multiple mappings exist between the query and reference genome.
For **secondary** alignments, these are overlapping and might indicate **duplications**. For **supplementary** alignments they are instead usually non-overlapping and might indicate **inversions/genome rearrangements**.

Figures for these alignments are stored in the `figures/vial_XX/non_primary` folder.
- `secondary_t_YY.pdf`: similar to a dotplot, where lines represent regions of the genome that have overlapping mapping. Black/red lines correspond to same-direction/inverted reads.
- `supplementary_t_YY.pdf`: scatter-plot, the x-y coordinates of each dot represent a link between two regions of the same query that map to region x and y of the reference in a chimeric read. Black/red dots correspond to same-direction/inverted reads. Only pairs of alignments that are non-overlapping on the query are selected.
- `secondary_vs_t.pdf` and `supplementary_vs_t.pdf`: histogram of starting positions of supplementary and secondary reads for the different timepoints.