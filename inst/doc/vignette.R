## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  library(devtools)
#  devtools::install_github("https://github.com/bowmanr/scDNA/tree/BB_dev",force=TRUE)

## ----message=FALSE,echo=FALSE,eval=FALSE--------------------------------------
#  library(scDNA)
#  library(googledrive)
#  library(dplyr)
#  library(Homo.sapiens)
#  library(VariantAnnotation)
#  library(GenomicRanges)
#  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  library(ggplot2)
#  

## ----load in file,eval=FALSE--------------------------------------------------
#  sample_file<- "~/CodingCamp/Projects/scDNA_demo/Sample1962.dna+protein.h5"
#  #sample_file<-c("/Users/bowmanrl/Desktop/Review/MSK24.cells.loom")

## ----variant identification,eval=FALSE----------------------------------------
#  #need to install SpareMatrixStats
#  variant_output<-variant_ID(file=sample_file,
#                             txdb="MSK_RL",
#                              GT_cutoff=0,
#                               VAF_cutoff=0)
#  
#  
#  genes_of_interest <- c("IDH2","NRAS","NPM1","TET2","FLT3","IDH1")

## ----variant selection orders,fig.align='center',fig.height=5,fig.width=5,eval=FALSE----
#  variants_of_interest<-variant_output%>%
#                          distinct()%>%
#                            dplyr::filter(Class=="Exon")%>%
#                            dplyr::filter(VAF>0.01)%>%
#                            dplyr::filter(genotyping_rate>85)%>%
#                            dplyr::filter(!is.na(CONSEQUENCE)&CONSEQUENCE!="synonymous")%>%
#                            dplyr::filter(SYMBOL%in%genes_of_interest)%>%
#                            dplyr::arrange(desc(VAF))%>%
#                            dplyr::slice(c(1,2,3))
#  

## ----sce construction,cache=TRUE,eval=FALSE-----------------------------------
#  sce<-tapestri_h5_to_sce(file=sample_file,
#                        variant_set = variants_of_interest,
#                        GT_cutoff=90,
#                        VAF_cutoff=0.01,
#                        DP_cutoff=10,
#                        GQ_cutoff=20,
#                        AF_cutoff=20)

## ----eval=FALSE---------------------------------------------------------------
#  sce<-enumerate_clones(sce, replicates = 500)
#  sce<-compute_clone_statistics(sce)

## ----clonograph,fig.align='center',fig.height=5,fig.width=5,eval=FALSE--------
#  clonograph(sce)

## ----subset clonograph,fig.align='center',fig.height=5,fig.width=5,eval=FALSE----
#  sce_subset<-select_clones(sce,
#                            ADO_cut=0.10,
#                            GQ_cut=30,
#                            DP_cut=10,
#                            select_exact=FALSE)
#  clonograph(sce,complete_only=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(compiler)
#  library(foreach)
#  library(doParallel)
#  library(visNetwork)
#  library(igraph)
#  #enableJIT(3)
#  sce<-trajectory_analysis(sce)
#  start_state ="0_0_0" # Change based on problem and number of variants we are using
#  goal_state = "1_1_1" # Change based on problem and number of varianrs we are using``

## ----eval=FALSE---------------------------------------------------------------
#  sce<-get_own_path(sce,start_state,goal_state)
#  visualize_any_optimal_path(sce,start_state,goal_state)
#  
#  print(sce@metadata$Trajectories$WT_to_last_possible_clone)
#  

## ----Visualize Whole Network,eval=FALSE---------------------------------------
#  vis_full_net<-visualize_full_network(sce)
#  vis_full_net

## ----visualize most likely order,fig.align='center',fig.height=5,fig.width=5,eval=FALSE----
#  vis_best_WT_to_dominant_clone<-visualize_WT_dominant_clone(sce)
#  vis_best_WT_to_dominant_clone
#  print(sce@metadata$Trajectories$WT_to_dominant_clone)

## ----visual all possible orders,fig.align='center',fig.height=5,fig.width=5,eval=FALSE----
#  vis_all_WT_to_dominant_clone<-visualize_all_WT_dominant_clone(sce)
#  vis_all_WT_to_dominant_clone
#  print(sce@metadata$Trajectories$All_WT_to_dominant_clone)
#  

## ----match clonal graph from WT,fig.align='center',fig.height=5,fig.width=5,eval=FALSE----
#  given_state <-"0"
#  match_clonal_graph(sce,given_state)
#  print(sce@metadata$Trajectories$WT_to_Observed_paths)
#  

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  library(dsb)
#  library(Seurat)
#  library(SingleCellExperiment)
#  droplet_metadata<- extract_droplet_size(sample_file,sce=sce_subset)
#  
#  # Try changing bins=1000 to see overlap empty droplets and real ones
#  ggplot(droplet_metadata, aes(x = dna_size, y = protein_size,color=Droplet_type )) +
#                    theme_bw() +
#                    geom_density_2d(bins=1000)+
#                    geom_hline(yintercept =c(1.5,5),lty=2)+
#                    geom_vline(xintercept =c(0.1,1.5),lty=2)

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  background_droplets<-droplet_metadata%>%
#                            dplyr::filter(Droplet_type=="Empty")%>%
#                            dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
#                            pull(Cell)
#  
#  sce<-normalize_protein_data(file=sample_file,
#                               metadata=droplet_metadata,
#                               sce=sce_subset,
#                               method=c("dsb","CLR"),
#                               detect_IgG=TRUE,
#                               background_droplets=background_droplets)
#  
#  altExp(sce)

## ----eval=FALSE---------------------------------------------------------------
#  library(Biobase)
#  library(flowCore)
#  library(flowViz)
#  names(assays(altExp(sce)))
#  fcs_export(sce,
#             slot="CLR_norm",
#             save_path="/Users/bowmanrl/Desktop/sample_CLR.fcs")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  library(Seurat)
#  
#  s<-Seurat::as.Seurat(x = altExp(sce),counts="Protein",data="CLR_norm")
#  s<-SeuratObject::RenameAssays(object = s, originalexp = 'Protein')
#  
#  colnames(s@meta.data)
#  colnames(s@meta.data)[12]<-"TET2"
#  colnames(s@meta.data)[13]<-"NPM1"
#  colnames(s@meta.data)[14]<-"NRAS"
#  colnames(s@meta.data)

## ----eval=FALSE---------------------------------------------------------------
#  s <- Seurat::ScaleData(s, assay = "Protein")
#  s <- Seurat::FindVariableFeatures(s,assay = "Protein")
#  s <- Seurat::RunPCA(s, features = VariableFeatures(object = s))
#  s <- Seurat::FindNeighbors(s, dims = 1:10)
#  s <- Seurat::FindClusters(s, resolution = 0.5)
#  # cluster and run umap (based directly on dsb normalized values without isotype controls)
#  
#  
#  s <- FindNeighbors(object = s,
#                     dims = NULL,
#                     k.param = 30)
#  
#  # direct graph clustering
#  s <- FindClusters(object = s,
#                    assay = "Protein",
#                    resolution = 1,
#                    algorithm = 3,
#                    verbose = TRUE)
#  
#  # umap for visualization only; (this is optional)
#  s = RunUMAP(object = s,
#              seed.use = 1990,
#              min.dist = 0.2,
#              n.neighbors = 50,
#              features=s@assays$Protein@var.features,
#              verbose = TRUE)
#  
#  DimPlot(s,reduction = "umap",group.by ="Clone" )
#  
#  FeaturePlot(s,features = c("CD19"))
#  Idents(s)<-"Clone"
#  Seurat::FindMarkers(s,ident.1 = "0_0_0",ident.2 = "1_1_1")
#  
#  Seurat::RidgePlot(s,features = "CD3")

