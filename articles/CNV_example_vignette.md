# CopyNumberVariation_tutorial

## CNV Tutorial for scDNA

In this tutorial we will go through an example of how we identify the
CNV properties of samples to perform analysis.

``` r
library(scDNA)
library(dplyr)
library(Homo.sapiens)
library(VariantAnnotation)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(SingleCellExperiment)
library(Seurat)
library(dsb) 
```

We will look into a TP53 sample with a identified TP53.Q331\* mutations
and complex karyotyping of del5q and a gain of chr8. First we start with
our standard scDNA pipeline.

``` r
sample_file = "../H5/CNV_example.dna+protein.h5"

variant_output<-variant_ID(file=sample_file,
                           panel="hg_38",
                            GT_cutoff=0,
                             VAF_cutoff=0)
genes_of_interest <- c("TP53")

variants_of_interest<-variant_output%>%
                        distinct()%>%
                          dplyr::filter(Class=="Exon")%>%
                          dplyr::filter(VAF>0.01)%>%
                         # dplyr::filter(genotyping_rate>85)%>%
                          dplyr::filter(!is.na(CONSEQUENCE)&CONSEQUENCE!="synonymous")%>% 
                          dplyr::filter(SYMBOL%in%genes_of_interest)%>%   
                          dplyr::arrange(desc(VAF))%>%
dplyr::group_by(final_annot)%>%
                          dplyr::mutate(final_annot = if(n( ) > 1) {paste0(final_annot,"_",LETTERS[row_number()]) } 
                             else {paste0(final_annot)})%>%
   ungroup()%>%
  dplyr::slice(c(1))

sce<-tapestri_h5_to_sce(file=sample_file,
                      variant_set = variants_of_interest,
                      GT_cutoff=90, 
                      VAF_cutoff=0.01,
                      DP_cutoff=10,
                      GQ_cutoff=20,
                      AF_cutoff=20)
sce<-enumerate_clones(sce, replicates = 500)
sce<-compute_clone_statistics(sce,skip_ploidy=FALSE)
```

Afterwards we want to run the standard pipeline to go from an sce object
to a seurat object and follow standard pipelines to identify cell types
in seurat

``` r
droplet_metadata<- extract_droplet_size(sce)

background_droplets<-droplet_metadata%>%
  dplyr::filter(Droplet_type=="Empty")%>%
  dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
  pull(Cell)

sce_subset<-normalize_protein_data(sce=sce,
                                   metadata=droplet_metadata,
                                   method=c("dsb","CLR"),
                                   detect_IgG=TRUE,
                                   background_droplets=background_droplets)

altExp(sce_subset)
s<-Seurat::as.Seurat(x = altExp(sce_subset),counts="Protein",data="CLR_norm")
s<-SeuratObject::RenameAssays(object = s, originalexp = 'Protein')

s <- Seurat::ScaleData(s, assay = "Protein")
s <- Seurat::FindVariableFeatures(s,assay = "Protein")
s <- Seurat::RunPCA(s, features = VariableFeatures(object = s))

s <- FindNeighbors(object = s,
                   dims = NULL,  
                   k.param = 20,
                   annoy.metric = "cosine")
# direct graph clustering 
s <- FindClusters(object = s, 
                  assay = "Protein",
                  resolution = 0.5, # smaller number means less clusters and bigger number means more clusters
                  algorithm = 3, 
                  verbose = TRUE)
s = RunUMAP(object = s, 
            seed.use = 1990, 
            min.dist = 0.2, 
            n.neighbors = 50, 
            features=s@assays$Protein@var.features,
            verbose = TRUE)


VlnPlot(s,idents = c("0","1","2","3","5","6"),features = c("CD34","CD117","CD123","CD11c","CD303","CD3"),pt.size = 0.001)
```

Differential marker expression to identify the cell type clusters.

``` r
marker_check = FindAllMarkers(s)

marker_check%>%group_by(cluster)%>%slice_head(n = 5)
```

Cell type labeling from Seurat protein content. This is currently
determined manually.

``` r
s<-RenameIdents(s,`0`="HSPC_CDP",`1`="NK",`2`="HSPC_MPP",`3`="Pre-DC",`4`="CD8+ T",`5`="Doublets",`6`="Granulocytes",`7`="CD4+ T")

s@meta.data$cell_type <-s@active.ident
s@meta.data$cell_type <-factor(s@meta.data$cell_type,levels=c("HSPC_MPP","HSPC_CDP","Pre-DC","Granulocytes","Doublets","NK","CD8+ T","CD4+ T"))
```

Let’s now get a dataframe that consists of cell type with their barcode.
When getting our CNV we will call a function called readDNA_CN_H5(). The
inputs for this function include sce object and reference cell barcodes.
We will use the T cell groups as our reference cells like below due to
their WT. the default for this setting is to use all cells. Note, when
running the compute_clone_statistics() we call readDNA_CN_H5() with
default settings inside.

``` r
meta_data_data_frame_cell_type <- FetchData(s,vars = "cell_type")

sce=readDNA_CN_H5(sce = sce,
                  reference_cells = meta_data_data_frame_cell_type%>% 
                    dplyr::filter(cell_type=="CD8+ T"|cell_type=="CD4+ T")%>%
                    rownames)
```

Now we can look into differences for ploidy. (See tutorial about cell
type VAFs if interested). We access the CNV data through the alternate
Experiment, altExp() function. We show how to extract the ploidy and
adapt it into a dataframe for further analysis.

``` r
CNV_sce = altExp(sce,e="CNV")
CNV_df<-CNV_sce@metadata$full_ploidy%>%
  dplyr::filter(!is.na(ploidy)&ploidy!=-Inf&ploidy!=Inf)%>%
  dplyr::arrange(CHROM,end_pos)%>%
  dplyr::select(amplicon,barcode,depth,CHROM,end_pos,cell_subset,ploidy,norm_frac,fraction)%>%
  dplyr::inner_join(sce@metadata$NGT_with_missing,by=c("barcode"="Cell"))%>%
  dplyr::group_by(amplicon)%>%
  dplyr::mutate(med_ploidy=median(ploidy))%>%
  dplyr::ungroup()%>%
  inner_join(meta_data_data_frame_cell_type%>%tibble::rownames_to_column(var="barcode"),by="barcode")%>%
  dplyr::mutate(cell_type=factor(cell_type,levels=c("HSPC_MPP","HSPC_CDP","Pre-DC","Granulocytes","Doublets","NK","CD8+ T","CD4+ T")))
```

Let’s plot the ploidy to compare chromosomes along cell type and
mutation status. We will see myeloid cells and TP53-mutation status
alter the complex karyotyping.

``` r
#Chromosome and cell type
gg_chrom_vs_celltype<-ggplot(CNV_df%>%filter(CHROM%in%c(5,8,13,4)),aes(x=CHROM,y=(2^(1+ploidy)),fill=cell_type))+
  geom_boxplot()+geom_hline(yintercept = 2,linetype='dashed')+
  coord_flip()+
  scale_fill_manual(values=pals::kelly(n=9)[-1])+
  xlab("Chromosome")+
  ylab("Number of Measured Chromosomes")+
  ylim(c(0,4))+
  ggpubr::theme_classic2(base_size=12)+
  theme()

gg_chrom_vs_mut<-ggplot(CNV_df%>%filter(CHROM%in%c(5,8,13,4)),aes(x=CHROM,y=(2^(1+ploidy)),fill=Clone))+
  geom_boxplot()+geom_hline(yintercept = 2,linetype='dashed')+
  #geom_jitter(size=0.2,alpha=0.2)+
  coord_flip()+
  scale_fill_manual(values=c("0"=RColorBrewer::brewer.pal(7,"Reds")[1],"1"=RColorBrewer::brewer.pal(7,"Reds")[3], "2"=RColorBrewer::brewer.pal(7,"Reds")[6],"3"="grey50"),label =c("WT","HET","HOM","DEL"))+
  xlab("Chromosome")+
  ylab("Number of Measured Chromosomes")+
  ylim(c(0,4))+
  ggpubr::theme_classic2(base_size=12)+
  guides(fill=guide_legend(title="TP53 Q331*"))

cowplot::plot_grid(plotlist = list(gg_chrom_vs_celltype,gg_chrom_vs_mut),align = "h",rel_widths = c(0.5,0.5))
```
