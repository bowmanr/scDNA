# The main function builds an adjacency list of the theoretical order of mutations. The same MDP mask can be reused different mutations

This contains 1 main function and 3 helper functions. The 3 helper
functions aim to:

1.  turn our states to Ternary variables (can be 0(WT), 1(Het), or
    2(Hom))

2.  a function that finds the difference in states that are \< 1
    mutation away

3.  a function that does bit logic subtraction for the ternary variables

## Usage

``` r
BuildMDP(num_mutations, use_ADO = FALSE)
```

## Arguments

- num_mutations:

  The number of variants we are using (automatically obtained from the
  Architecture)
