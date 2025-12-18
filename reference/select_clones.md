# Select clones of interest on the basis QC metrics

Select clones of interest on the basis QC metrics

## Usage

``` r
select_clones(
  sce,
  ADO_cut = 0.1,
  GQ_cut = 30,
  DP_cut = 10,
  select_exact = FALSE
)
```

## Arguments

- sce:

  a SingleCellExperiment object containing the Clones,Architecture, and
  NGT matrix

- ADO_cut:

  maximum median allele dropout

- GQ_cut:

  mimimum median gene quality score

- DP_cut:

  minimum median sequencing depth

- select_exact:

  default=FALSE, can supply a vector of clones of interest to bypass QC
  filtering

## Value

a subset of the sce object
