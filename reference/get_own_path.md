# A way to navigate the MDP for any given starting root node to leaf node.

A way to navigate the MDP for any given starting root node to leaf node.

## Usage

``` r
get_own_path(sce, start_name, goal_name)
```

## Arguments

- sce:

  a singleCellExperiment object

- start_name:

  a clone_code of starting point. make sure it is less than goal_name

- goal_name:

  a clone_code for the end point to get to. make sure it is larger than
  goal_name

## Value

It augments the sce object in the Trajectories list in the meta data by
adding another path
