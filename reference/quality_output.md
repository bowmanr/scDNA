# Produce long form quality metrics for each cell-variant pair for read depth, allele frequency, and genotype quality

Produce long form quality metrics for each cell-variant pair for read
depth, allele frequency, and genotype quality

## Usage

``` r
quality_output(
  file,
  filter = TRUE,
  input_variants,
  input_cells,
  NGT,
  DP_cut = 10,
  AF_cut = 20,
  GQ_cut = 20
)
```

## Arguments

- file:

  path to H5 file

- filter:

  default=TRUE, determines whether function should output filtered or
  unfiltered long form matrix containing per cell-variant pair of NGT,
  DP, AF, GQ metrics

- input_variants:

  variants of interest to collect QC data

- input_cells:

  cells of interest to subset on

- NGT:

  NGT matrix

- DP_cut:

  read depth cutoff

- AF_cut:

  allele frequency cutoff

- GQ_cut:

  genotype quality cutoff

## Value

filtered NGT matrix containing NAs for cells that do not pass QC
