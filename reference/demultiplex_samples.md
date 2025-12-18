# Demultiplex wrapper function to handle splitting samples

Demultiplex wrapper function to handle splitting samples

## Usage

``` r
demultiplex_samples(
  sce,
  sensitivity_threshold = c(0.01, 1e-04),
  expected_samples = 5
)
```

## Arguments

- sce:

  a SingleCellExperiment object containing the Clones,Architecture, and
  NGT matrix

- sensitivity_threshold:

  vector of sensitivity thresholds for rows (variants) and columns
  (cells)

- expected_samples:

  the number of clusters or samples you expect

## Value

a dataframe of the cells assigned to a specific cluster.
