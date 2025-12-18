# BSCITE Style single trajectory Plot

BSCITE Style single trajectory Plot

## Usage

``` r
trajectory_of_interest_BSCITE_format(
  sce,
  trajectory = sce@metadata$Trajectories[[1]],
  save_filename = NULL
)
```

## Arguments

- sce:

  SingleCellExperiment object containing NGT matrix for clone
  identification.

- trajectory:

  Put in entire Trajectory you are interested in such as
  'sce(at)metadata\$Trajectories\[1\]'

- save_filename:

  file location to save file. Currently, must end it with .html
  extension.
