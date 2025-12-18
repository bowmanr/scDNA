# Enumerate clones

Enumerate clones

## Usage

``` r
enumerate_clones(sce, replicates = 100)
```

## Arguments

- sce:

  SingleCellExperiment object containing NGT matrix for clone
  identification.

- replicates:

  number of bootstrapping replicates

## Value

updated sce containing a table of clones, modified NGT matrix, and
clonal architecture.
