library(devtools)
library(tidyverse)
library(fs)

###create_package("~/Projects/R_packages/scDNA")
use_git()
use_mit_license("Bobby Bowman")

setwd("~/Projects/R_packages/scDNA")

use_r("read_tapestri_h5_NGT")
use_r("quality_filter_NGT")
use_r("quality_output")
use_r("annotate_variants")
use_r("TI_mimic_filter")
use_r("active_NGT_filter")
use_r("generate_lowQC_matrix")
use_r("enumerate_clones")
use_r("clonograph")
use_r("clone_QC")
use_r("select_clones")
use_r("read_tapestri_h5_protein")
use_r("extract_droplet_size")
use_r("normalize_protein_data")
use_readme_rmd()

use_devtools()

load_all()

document()
check()


use_package("TxDb.Hsapiens.UCSC.hg19.knownGene")
setwd("/Users/bowmanr/Projects/scDNA/")
file<-("./data/Sample1962.dna+protein.h5")

NGT<-read_tapestri_h5_NGT("/Users/bowmanr/Projects/scDNA/data/Sample1962.dna+protein.h5")


filtered_NGT<-quality_filter_NGT(file="/Users/bowmanr/Projects/scDNA/data/Sample1962.dna+protein.h5",
                                 NGT=NGT,
                                 DP_cut=10,
                                 AF_cut=20,
                                 GQ_cut=20)

annotation_key <-read.csv("/Users/bowmanr/Projects/scDNA/data/annotation_key.csv")
hg19refseq_txdb<-loadDb(file="/Users/bowmanr/Projects/scDNA/data/hg19refseq_txdb.sqlite")
banned <-read.csv("/Users/bowmanr/Projects/scDNA/data/banned_list.csv")

annotation_key%<>%inner_join(select(hg19refseq_txdb,
                                    keys=annotation_key$ccds_id,
                                    columns=c("TXID","TXNAME"),
                                    keytype = "TXNAME"),
                             by=c("ccds_id"="TXNAME"))%>%
                  mutate(TXID=as.character(TXID))

final_mutation_info<- annotate_variants(file="/Users/bowmanr/Projects/scDNA/data/Sample1962.dna+protein.h5",
                             annotation_key=annotation_key,
                             txdb=hg19refseq_txdb,
                             banned=banned,#[,1],
                             NGT=filtered_NGT)
