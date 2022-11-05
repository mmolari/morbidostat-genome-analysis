# non-primary mappings

Othe than primary mappings, a `sam`/`bam` file can also contain:
- unmapped reads
- secondary alignments
- supplementary alignments

Example for one particular query that has both secondary and supplementary aligments:
| map-n | read-id                              | flag | fwd   | ref_len | qry_len | ref-start | ref-end | qry-start | qry-end |
| ----- | ------------------------------------ | ---- | ----- | ------- | ------- | --------- | ------- | --------- | ------- |
| 10461 | 00813e02-566a-405a-8adc-4ff12225e8c6 | 16   | False | 7391    | 8822    | 326310    | 333701  | 0         | 7357    |
| 38691 | 00813e02-566a-405a-8adc-4ff12225e8c6 | 256  | True  | 1328    | 0       | 1368099   | 1369427 | 145       | 1467    |
| 80164 | 00813e02-566a-405a-8adc-4ff12225e8c6 | 256  | True  | 1329    | 0       | 3275297   | 3276626 | 145       | 1467    |
| 95819 | 00813e02-566a-405a-8adc-4ff12225e8c6 | 2048 | True  | 1342    | 1336    | 3867913   | 3869255 | 0         | 1336    |

## unmapped reads

These are reads that have not been mapped to any part of the reference. Our pipeline stores a list in the `pileup/unmapped.csv` file, with entries:
| read                                 | len | avg. qscore        | flag |
| ------------------------------------ | --- | ------------------ | ---- |
| ffe7d278-7558-4fab-bdfe-ac8a646aca93 | 306 | 18.202614379084967 | 4    |
| bb873e5c-36e7-49dc-a442-b2cf9fdf6892 | 289 | 48.20415224913495  | 4    |
| e6172fcf-d15a-4f41-bdf0-7d6a014d7811 | 277 | 22.938628158844764 | 4    |
| ...                                  | ... | ...                | ...  |

Moreover the unmapped reads are extracted and saved in `pileup/unmapped.fastq.gz`.

## non-primary alignments

### secondary alignments

These corresponds to reads for which a particular sub-region is mapped in multiple positions (i.e. **multiple alignments / duplicated regions**). For these reads the *sequence* and *quality* is not reported in the SAM file (?).

These can have flags:
- 256 `0b100000000`, secondary
- 272 `0b100010000`, secondary + reverse-complemented

### supplementary alignments

These corresponds to reads for which different parts of the query map to different parts of the reference (i.e. **chimeric alignments/genome rearrangements**).

These can have flags:
- 2048 `0b100000000000`, supplementary
- 2064 `0b100000010000`, supplementary + reverse-complemented

### analysis

In our pipeline information on secondary and supplementary alignments is stored in `pileup/non_primary.csv`. This is a dataframe with summary information on all reads that have secondary or supplementary alignments (also includes primary alignments).