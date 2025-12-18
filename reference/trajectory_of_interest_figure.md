# Single Trajectory visualization of interest.

Single Trajectory visualization of interest.

## Usage

``` r
trajectory_of_interest_figure(
  sce,
  trajectory = sce@metadata$Trajectories[[1]],
  save_filename = NULL,
  tree_flag = FALSE
)
```

## Arguments

- sce:

  SingleCellExperiment object containing NGT matrix for clone
  identification.

- trajectory:

  This is any trajectory you wanted to create in the
  metadata\$Trajectory

- save_filename:

  file location to save file. Currently, must end it with .html
  extension.

- tree_flag:

  set TRUE if you want a tree layout, default is FALSE and freeform to
  move around
