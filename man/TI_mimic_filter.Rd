% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TI_mimic_filter.R
\name{TI_mimic_filter}
\alias{TI_mimic_filter}
\title{Function to mimic basic variant and cell cutoffs}
\usage{
TI_mimic_filter(
  NGT_filter,
  variants_of_interest,
  required_variants,
  required_cells,
  variant_score_cutoff = 50,
  cell_score_cutoff = 50
)
}
\arguments{
\item{NGT_filter}{required, filtered NGT matrix containing 0,1,2 and NA for cells on rows, and variants in columnns}

\item{variants_of_interest}{required, universe of potential variants that have been prefiltered}

\item{required_variants}{optional, vector of variants that cannot be removed}

\item{required_cells}{optional, vector of cells that cannot be removed}

\item{variant_score_cutoff}{required, convergence value for fraction of cells that a given variant is successfully genotyped for scale is 0-1}

\item{cell_score_cutoff}{required,  convergence value for fraction of variants that a given cell is successfully genotyped for scale is 0-1}
}
\value{
NGT matrix
}
\description{
Function to mimic basic variant and cell cutoffs
}
