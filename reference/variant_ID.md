# Variant identification and frequency tallies The DNA variants for each cell are pulled from the specified H5 files. Each variant for each cell is genotyped to be WildType, Heterozygous, Homozygous or Missing. The genotyping rate is determined by taking WT+Het+Hom over total cell calls (including missing). The VAF is determined by number of allele copies we see in a weighted sum. A filter is applied to both of these calculatins to include or exclude variants of interest. These are then annotated to include the variant information such as gene name, nucleotide location, and short amino acid changes.

Variant identification and frequency tallies The DNA variants for each
cell are pulled from the specified H5 files. Each variant for each cell
is genotyped to be WildType, Heterozygous, Homozygous or Missing. The
genotyping rate is determined by taking WT+Het+Hom over total cell calls
(including missing). The VAF is determined by number of allele copies we
see in a weighted sum. A filter is applied to both of these calculatins
to include or exclude variants of interest. These are then annotated to
include the variant information such as gene name, nucleotide location,
and short amino acid changes.

## Usage

``` r
variant_ID(
  file,
  panel = NULL,
  GT_cutoff = 0,
  VAF_cutoff = 0,
  demultiplexed = NULL
)
```

## Arguments

- file:

  path to the h5 file

- panel:

  name of prebuilt panel/txdb

- GT_cutoff:

  Fraction of cells that are successfully genotyped for initial
  filtering (default 0.2, meaning 20%)

- VAF_cutoff:

  Fraction of cells that are mutated for initial filtering of variants
  (default 0.005, meaning 0.05%)

- demultiplexed:

  cell to cluster dataframe, often is set to NULL

## Value

A dataframe with each variant on a row, and tally of the number of cells
that are WT, Het, Hom or missing for a mutation. Calculated VAF and
gentoyping frequency is also provided. If multiple samples are present
in the h5 file, a list object will be returned with each sample as an
entry in the list
