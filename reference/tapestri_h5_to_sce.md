# Import Tapestri H5 data and extract genotype matrix

Import Tapestri H5 data and extract genotype matrix

## Usage

``` r
tapestri_h5_to_sce(
  file,
  GT_cutoff = 35,
  VAF_cutoff = 5,
  DP_cutoff = 10,
  GQ_cutoff = 30,
  AF_cutoff = 25,
  variant_set = NULL,
  demultiplex_cells = NULL,
  protein = TRUE
)
```

## Arguments

- file:

  path to the h5 file

- GT_cutoff:

  Fraction of cells that are successfully genotyped for initial
  filtering (default 0.2, meaning 20%)

- VAF_cutoff:

  Fraction of cells that are mutated for initial filtering of variants
  (default 0.005, meaning 0.05%)

- DP_cutoff:

  minimum number of reads necessary for a reliable genotype call in a
  single cell (default: 10)

- GQ_cutoff:

  minimum genotype quality necessary for a reliable genotype call in a
  single cell (default: 30)

- AF_cutoff:

  Deviation from 0, 50, or 100% for a reliable call of WT, Het or Hom
  respectively (default: 25)

- variant_set:

  character vector of variants to be included in the format output by
  mission bio.

- demultiplex_cells:

  This parameter is for demultiplexing cell names to use, often set to
  NULL

- protein:

  logical, whether protein data should be included. default=TRUE

- return_variants_only:

  logical,

## Value

a single cell experiment object containing the genotyping matrix, allele
frequency table, annotation table,
