# possible improvements

- [ ] Rename records in original fasta/genbank assembled genome as `vialXX_timeYY_clZZ`?
- [x] Group together correlated deletion trajectories that span adjacen intervals.
- [ ] parallelize insertion fisher test evaluation.
- [ ] because of quality filtering the number of insertions might be greater than the number of reads in the pileup. To solve for this in the insertion analysis we decrease the number of insertions in this case, to have an insertion frequency not greater than one. However this might bias the estimation for insertion frequency on the rest of the sites. A proper solution would require also keeping information on the unfiltered number of reads (also an unfiltered pileup).
- [x] Select trajectories based on delta (max - min) frequency on timepoints with high confidence.
- [x] Use secondary/supplementary reads to find duplicated/chimeric region bridges.
- [ ] Make the pipeline less reliant on folder structure. Pass the files directly in channels.
- [ ] Merge "pileup" and "extract" workflows.
- [ ] Move secondary/supplementary alignment plots from pileup to plot workflow. Maybe add extraction script to create a single database?
- [ ] Update documentation.