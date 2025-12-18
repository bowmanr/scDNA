# This function generates the Copy Number by determining the ploidy of each mutation

This function generates the Copy Number by determining the ploidy of
each mutation

## Usage

``` r
readDNA_CN_H5(sce, reference_cells = NULL)
```

## Arguments

- sce:

  the single cell experiment object

- reference_cells:

  set cells that will be used to determine allele dropout and false
  postive rates
