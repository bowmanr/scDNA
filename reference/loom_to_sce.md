# Title

Title

## Usage

``` r
loom_to_sce(
  file,
  GT_cutoff = 35,
  VAF_cutoff = 5,
  DP_cutoff = 10,
  GQ_cutoff = 30,
  AF_cutoff = 25,
  variant_set = NULL
)
```

## Arguments

- file:

  path to h5 file

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

## Value

list
