% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_clones.R
\name{select_clones}
\alias{select_clones}
\title{Select clones of interest on the basis QC metrics}
\usage{
select_clones(
  final_sample_summary = final_sample_summary,
  ADO_cut = 0.1,
  GQ_cut = 30,
  DP_cut = 10,
  select_exact = FALSE
)
}
\arguments{
\item{ADO_cut}{maximum median allele dropout}

\item{GQ_cut}{mimimum median gene quality score}

\item{DP_cut}{minimum median sequencing depth}

\item{select_exact}{default=FALSE, can supply a vector of clones of interest to bypass QC filtering}
}
\value{
final_sample_summary
}
\description{
Select clones of interest on the basis QC metrics
}
