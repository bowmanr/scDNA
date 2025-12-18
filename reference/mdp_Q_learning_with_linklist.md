# This file run Reinforcement Learning (model-free Q-learning) to evaluate the most likely mutation paths from the data

This file run Reinforcement Learning (model-free Q-learning) to evaluate
the most likely mutation paths from the data

## Usage

``` r
mdp_Q_learning_with_linklist(df, discount, N)
```

## Arguments

- df:

  this is the adjacency link link to

- discount:

  the factor determines how much future mutations impact current
  actions. 0 only looks at next step, close to 1 will have a very long
  horizon time

- N:

  Number of iterations, also called episodes to explore the MDP

## Value

A dataframe that contains the MDP adjacency list, normalized Q-values
and mutations taken
