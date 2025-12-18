# Annotate variants of interest This function takes in h5 files to extract DNA information through given TXDB files (primarily hg19) The gene names are mapped for their specific nucleotide positions. The starting and ending positions for the chromosome are listed. The consequence of the mutation such as synonymous, nonsynonymous, as well as if this is a coding region, on the exon boundry or intronic. Short amino acid changes are also labeled. The variance matrix also extracts the amplicon that the variant is found on.

Annotate variants of interest This function takes in h5 files to extract
DNA information through given TXDB files (primarily hg19) The gene names
are mapped for their specific nucleotide positions. The starting and
ending positions for the chromosome are listed. The consequence of the
mutation such as synonymous, nonsynonymous, as well as if this is a
coding region, on the exon boundry or intronic. Short amino acid changes
are also labeled. The variance matrix also extracts the amplicon that
the variant is found on.

## Usage

``` r
annotate_variants(file, panel = NULL, select_variants = NULL)
```

## Arguments

- file:

  path to h5 file to pull out relevant DNA information

- select_variants:

  variants of interest, default is all variants

## Value

variant annotation matrix
