# Plotting clonographs

Plotting clonographs

## Usage

``` r
clonograph(sce, complete_only = FALSE, color_pal = "Reds", QC_stats = FALSE)
```

## Arguments

- sce:

  the singleCellExperiment Object with the clones, NGT, and
  architecture.

- complete_only:

  logical to decide whether to plot only complete cells or also low
  quality cells, default=FALSE

- color_pal:

  brewer.pal pallete selection.

- QC_stats:

  A true/false call to see other diagnostic plots under the clonograph
  (default=FALSE)

## Value

A clonograph figure
