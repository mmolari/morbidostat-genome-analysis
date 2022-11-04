## unmapped reads

## secondary alignment

Corresponds to **multiple alignments** (duplicated regions).

- sequence and quality not reported

Flags:
- 256 `0b100000000`, secondary
- 272 `0b100010000`, secondary + reverse-complemented


## supplementary alignment

Corresponds to **chimeric alignments** (genome rearrangements?)

Flags:
- 2048 `0b100000000000`, supplementary
- 2064 `0b100000010000`, secondary + reverse-complemented
