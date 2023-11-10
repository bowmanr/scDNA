
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDNA v1.0.1

<!-- badges: start -->
<!-- badges: end -->

The goal of scDNA R package is to provide a simple framework for
analyzing single cell DNA sequencing data. The current version primarily
focuses processing variant information on the Mission Bio Tapestri
platform. Functionality includes import of h5 files from Tapestri
pipeline, basic variant annotation, genotype extraction, clone
identification, and clonal trajectory inference. This package provides
wrappers for normalizing protein data for scDNA+Protein libraries for
downstream analysis.

## Installation

You can install the current version (1.0.1) of scDNA below

``` r
remotes::install_github("bowmanr/scDNA")
```

## Version Updates

### **v1.0.1**

- H5 files are now read using the rhdf5 package and stored into a
  [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
  container.

  - Merged h5 samples are identified and sample names are stored in
    colData(). Variant identification is ran separately and then merged.

  - Variant information is stored in rowData()

  - NGT matrix, clonal abundance, and clone architecture [familiar to
    previous versions](https://bowmanr.github.io/scDNA_myeloid/) can be
    found in the metadata.

- Variant identification and annotation is performed initially before
  reading in all the genotyping/QC data.

  - Transcript annotation matches [cannonical
    transcripts](https://docs.cbioportal.org/mutation-data-transcript-annotation/#transcript-assignment)
    used in the [cBio portal](https://www.cbioportal.org/).

  - To decrease variant location identification runtime, we created a
    custom TxDB object for the Clonal Evolution Panel from used
    [here](https://www.nature.com/articles/s41586-020-2864-x). If you
    have a different panel you can also use the [TxDB for hg19 from
    UCSC.](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html)
    Future versions will have local data for all panels from Mission
    Bio, as well as a simple script for generating a TxDB object for
    custom panels.

- Protein data is stored as an altExp() container within the container.

  - Wrappers for [DSB](https://github.com/niaid/dsb) and CLR
    normalization are provided. (CLR currently performed in
    [Seurat](https://satijalab.org/seurat/)).

  - Simple import into Seurat is demonstrated.

  - Export to FCS files with mutations and clone “completeness” provided
    as variables.

## Simple workflow

Identify all variants within a sample.

``` r
library(scDNA)
library(dplyr)
sample_file<- "test_file.h5"
variant_output<-variant_ID(file=sample_file,
                           panel="MSK_RL", # "UCSC" can be used for other panels
                           GT_cutoff=0,  # mimimum percent of cells where a successful genotyping call was made
                           VAF_cutoff=0) # mimimum variant allele frequency 
```

Identify mutations in genes of interest.

``` r
genes_of_interest <- c("IDH2","NRAS","NPM1","TET2","FLT3","IDH1")
variants_of_interest<-variant_output%>%
                          dplyr::filter(Class=="Exon")%>%
                          dplyr::filter(VAF>0.01)%>%
                          dplyr::filter(genotyping_rate>85)%>%
                          dplyr::filter(!is.na(CONSEQUENCE)&CONSEQUENCE!="synonymous")%>%
                          dplyr::filter(SYMBOL%in%genes_of_interest)%>%   
                          dplyr::arrange(desc(VAF))%>%
                          dplyr::slice(c(1:3)) # take the 3 most abundance mutations
```

Read in the data, enumerate clones, and compute statistics. Sample
statistics mirror that seen in Figure 1
[here](https://www.nature.com/articles/s41586-020-2864-x), and are
stored in the metadata.

``` r
sce<-tapestri_h5_to_sce(file=sample_file,variant_set = variants_of_interest)
sce<-enumerate_clones(sce)
sce<-compute_clone_statistics(sce)
```

Simple function for producing a graph in the style of Figure 1D from
[here](https://www.nature.com/articles/s41586-020-2864-x),

``` r
clonograph(sce)
```

<img src="images/Screen%20Shot%202023-09-19%20at%2010.07.16%20PM.png"
width="373" />

Function to perform Reinforcment Learning / MDP approach for clonal
trajectory as in Figure 3
[here](https://www.nature.com/articles/s41586-020-2864-x),

``` r
sce<-trajectory_analysis(sce)
```

Methods for protein normalization. Both dsb and CLR normalization can be
performed and stored in separate slots. We tend to have favor dsb so
far.

``` r
droplet_metadata<- extract_droplet_size(sce)
background_droplets<-droplet_metadata%>%
                          dplyr::filter(Droplet_type=="Empty")%>%
                          dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
                          pull(Cell)

sce<-normalize_protein_data(sce=sce,
                             metadata=droplet_metadata,
                             method=c("dsb","CLR"),
                             detect_IgG=TRUE,
                             background_droplets=background_droplets)
```

### Developments in progress:

1.  Cohort summarization
2.  Creating custom TxDB objects
3.  Improve trajectory plotting
4.  Sample demultiplexing (integration from [Robinson et
    al](https://www.biorxiv.org/content/10.1101/2022.09.20.508786v1.full),
    [github](https://github.com/RobinsonTroy/scMRD))
5.  Evaluation of allele dropout and clonal confidence
6.  Integration for copy number variation

### Ongoing investigation:

1.  Improving cell identification and distinction from empty droplets.
    1.  Doublet and dead cell identification
2.  Improve normalization for protein data.
    1.  Improve cell type identification based on immunophenotype
3.  Improvements to the MDP and RL.
