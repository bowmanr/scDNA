
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDNA

<!-- badges: start -->
<!-- badges: end -->

The goal of scDNA R package is for DNA+protein single cell DNA
sequencing analysis. This includes extracting data from Tapestri, to
produce genotype matrix (NGT), determine variant allele frequency(VAF),
summarize clone information, and determine the clonal trajectories.

## Installation

You can install the development version of scDNA like so:

``` r
devtools::install_github("bowmanr/scDNA")
```

    #> 
    #> ── R CMD build ─────────────────────────────────────────────────────────────────
    #> * checking for file ‘/private/var/folders/x3/kjrtpfnn0_gf6gqk86ls4kmm0000gp/T/RtmpJqBPJM/remotesc6be610ce9d7/bowmanr-scDNA-acd5013/DESCRIPTION’ ... OK
    #> * preparing ‘scDNA’:
    #> * checking DESCRIPTION meta-information ... OK
    #> * checking vignette meta-information ... OK
    #> * checking for LF line-endings in source and make files and shell scripts
    #> * checking for empty or unneeded directories
    #> * building ‘scDNA_0.1.0.9000.tar.gz’

Stay tuned for the following updates and things we are actively working
on next:

### Near term:

1.  FCS export
2.  Finalizing trajectory optimization
3.  Finalizing cohort characterization
4.  Creating your own TxDB tutorial

### Ongoing investigation:

1.  Improved cell identification.
2.  A new filtering and normalization for protein data.
3.  Diagnostics on your assays.
4.  Improvements to the MDP and RL:
5.  Introducing safe RL through risk-sensitivity
6.  Introducing probabilistic states and actions
7.  Improvements in performance by changing to a sparse representation
8.  Improved trajectory analysis and path finding. 1. longest pathway to
    mutation 2. free energy principle to find shortest path.
9.  improvements to MDP formalization for improved accuracy and
    efficiency.
10. Improvements to RL algorithm to speed up the learned model.
11. Updating presentation of figures.
