#' Import Tapestri H5 data and extract genotype matrix
#'
#' @param file path to the h5 file
#' @param GT_cutoff Fraction of cells that are successfully genotyped for initial filtering (default 0.2, meaning 20%)
#' @param VAF_cutoff Fraction of cells that are mutated for initial filtering of variants (default 0.005, meaning 0.05%)
#' @param DP_cutoff minimum number of reads necessary for a reliable genotype call in a single cell (default: 10)
#' @param GQ_cutoff minimum genotype quality necessary for a reliable genotype call in a single cell (default: 30)
#' @param AF_cutoff Deviation from 0, 50, or 100% for a reliable call of WT, Het or Hom respectively (default: 25)
#' @param variant_set character vector of variants to be included in the format output by mission bio.
#' @param protein logical, whether protein data should be included. default=TRUE
#' @param return_variants_only logical,
#'
#' @return a single cell experiment object containing the genotyping matrix, allele frequency table, annotation table,
#' @export
#'
#' @examples
tapestri_h5_to_sce<-function(file,
                    GT_cutoff=35,
                    VAF_cutoff=5,
                    DP_cutoff=10,
                    GQ_cutoff=30,
                    AF_cutoff=25,
                    variant_set=NULL,
                    final_mutation_set=NULL,
                    protein=TRUE,
                    return_variants_only=FALSE){
  
  
  if(is.null(variant_set)){
    print("No variants provided, running in discovery mode.
          If you want to take a look at variants only use the 'variant_ID' function, or set 'return_variants_only==TRUE'")
    total_variants <- scDNA::variant_ID(file=file,GT_cutoff=GT_cutoff,VAF_cutoff = VAF_cutoff)
    VAF_cut_variants<- total_variants%>%
      dplyr::pull(id)
    annotated_variants<- scDNA::annotate_variants(file,select_variants = NULL)
    total_variants[[1]]<-annotated_variants%>% 
      dplyr::inner_join(total_variants[[1]],by="id")%>%
      dplyr::filter(!is.na(SYMBOL))%>%
      dplyr::filter(!is.na(CONSEQUENCE)&CONSEQUENCE!="synonymous")%>%
      dplyr::arrange(desc(VAF))
    if(return_variants_only==TRUE){
      return(total_variants)
      break
    } 
  } 
  
  if(!is.null(variant_set)){
    # variant_set has all info, so we need to pull out the id for the h5 file
    VAF_cut_variants<- variant_set%>%
      dplyr::filter(genotyping_rate>=GT_cutoff)%>%
      dplyr::filter(VAF>=VAF_cutoff)%>%
      dplyr::arrange(desc(VAF))%>%
      pull(id)
    
    print(paste("Loading n=",length(VAF_cut_variants),"variants")) 
    print(paste("VAF Range:", variants_of_interest%>%reframe(range=range(VAF))%>%pull(range)%>%round(digits = 2)))
    
    print(paste("Input file:",file))
    sample_set<-rhdf5::h5read(file=file,name="/assays/dna_variants/metadata/sample_name")[1,]
    
    print(paste("Detected n =",length(sample_set),"sample(s):", paste(sample_set,sep=" ",collapse = " ")))
    
    print("Checking for Barcode Duplicates")
    all_barcodes<- rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
    dedup_barcodes<-all_barcodes[!(duplicated(all_barcodes) | duplicated(all_barcodes, fromLast = TRUE))]
    viable_barcodes<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%dedup_barcodes)
    print(paste(length(all_barcodes)-length(viable_barcodes), "duplicated barcodes detected and removed"))
  }
  
  VAF_cut_index<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")%in%VAF_cut_variants)

  print("Loading Allele Frequency Data")           
  AF_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/AF",index=list(VAF_cut_index,viable_barcodes))%>%
    `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
    data.frame()%>%
    dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
    tidyr::pivot_longer(cols=!id,values_to = "AF",names_to = "barcode")
  
  print("Loading Depth Data")
  print(paste("Depth Cutoff >",DP_cutoff))
  DP_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/DP",index=list(VAF_cut_index,viable_barcodes))%>%
    `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
    data.frame()%>%
    dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
    tidyr::pivot_longer(cols=!id,values_to = "DP",names_to = "barcode")
  
  
  print("Loading Genotype Quality Data")
  print(paste("Genotype quality cutoff >",GQ_cutoff))
  GQ_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/GQ",index=list(VAF_cut_index,viable_barcodes))%>%
    `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
    data.frame()%>%
    dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
    tidyr::pivot_longer(cols=!id,values_to = "GQ",names_to = "barcode")
  
  print("Loading Subsetted Genotype Information")
  NGT_data_subset<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(VAF_cut_index,viable_barcodes))%>%
    `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
    data.frame()%>%
    dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
    tidyr::pivot_longer(cols=!id,values_to = "NGT",names_to = "barcode")
  
  print("Final Filtering")
  sample_names<-data.frame("barcode"=rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)),    
                           "Sample"=rhdf5::h5read(file=file,name="/assays/dna_variants/ra/sample_name",index=list(viable_barcodes)))%>%
    dplyr::mutate(barcode=gsub("-","\\.",barcode))
  
  
  out<-list(AF_data,
            GQ_data,
            DP_data,
            NGT_data_subset)%>%
    purrr::reduce(dplyr::full_join, by=c("id","barcode"))%>%
    dplyr::full_join(sample_names,by="barcode")%>%
    dplyr::mutate(NGT_mask=ifelse(NGT==3,FALSE,TRUE))%>%
    dplyr::mutate(AF_mask=case_when( 
      (NGT==0 & AF<(50-AF_cutoff))~TRUE,
      (NGT==1 & (AF>(50-AF_cutoff)|AF<(50+AF_cutoff)))~TRUE,
      (NGT==2 & AF>(50+AF_cutoff))~TRUE,
      TRUE~FALSE
    ))%>%
    dplyr::mutate(DP_mask=ifelse(DP>DP_cutoff,TRUE,FALSE))%>%
    dplyr::mutate(GQ_mask=ifelse(GQ>GQ_cutoff,TRUE,FALSE))
  

  final_barcodes<-out$barcode%>%unique%>%{gsub("\\.","-",.)}
  final_barcodes_index<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%final_barcodes)
  
  print("Formating Single Cell Experiment")
  
  sce<-SingleCellExperiment::SingleCellExperiment(list(
    NGT=as.data.frame(tidyr::pivot_wider(out,
                           id_cols=id,
                           names_from=barcode,
                           values_fill=NA,
                           values_from = NGT))%>%
      dplyr::select(-id),
    NGT_mask=as.data.frame(tidyr::pivot_wider(out,
                                id_cols=id,
                                names_from=barcode,
                                values_fill=NA,
                                values_from = NGT_mask))%>%
      dplyr::select(-id),                              
    AF=as.data.frame(tidyr::pivot_wider(out,
                          id_cols=id,
                          names_from=barcode,
                          values_fill=NA,
                          values_from = AF))%>%
      dplyr::select(-id),
    AF_mask=as.data.frame(tidyr::pivot_wider(out,
                               id_cols=id,
                               names_from=barcode,
                               values_fill=NA,
                               values_from = AF_mask))%>%
      dplyr::select(-id),
    DP=as.data.frame(tidyr::pivot_wider(out,
                          id_cols=id,
                          names_from=barcode,
                          values_fill=NA,
                          values_from = DP))%>%
      dplyr::select(-id),
    DP_mask=as.data.frame(tidyr::pivot_wider(out,
                               id_cols=id,
                               names_from=barcode,
                               values_fill=NA,
                               values_from = DP_mask))%>%
      dplyr::select(-id),
    GQ=as.data.frame(tidyr::pivot_wider(out,
                          id_cols=id,
                          names_from=barcode,
                          values_fill=NA,
                          values_from = GQ))%>%
      dplyr::select(-id),
    GQ_mask=as.data.frame(tidyr::pivot_wider(out,
                               id_cols=id,
                               names_from=barcode,
                               values_fill=NA,
                               values_from = GQ_mask))%>%
      dplyr::select(-id)))
  # New addition of the annotation table to the single cell object.
  SummarizedExperiment::colData(sce)<-S4Vectors::DataFrame(sample=out%>%
                                                             dplyr::distinct(barcode,Sample)%>%
                                                             pull(Sample))
  rownames(sce)<-tidyr::pivot_wider(out,
                                    id_cols=id,
                                    names_from=barcode,
                                    values_fill=NA,
                                    values_from = NGT)%>%
    dplyr::pull(id)
  
  SummarizedExperiment::rowData(sce)<-S4Vectors::DataFrame(variant_set%>%
                                                             dplyr::rename(Widht=width,Strand=strand,Seqnames=seqnames,Start=start,End=end)%>%
                                                             dplyr::arrange(factor(id,levels=rownames(sce))))
  
 colnames(sce)<-tidyr::pivot_wider(out,
                                    id_cols=id,
                                    names_from=barcode,
                                    values_fill=NA,
                                    values_from = NGT)%>%
    dplyr::select(-id)%>%
    colnames()
  rownames(sce)<-rowData(sce)$final_annot

  print("Reordering sce rows (variants) based on bulk VAF")
  # this is newly added so we can get correct order for annotation later.
  sce<- sce[match(SummarizedExperiment::rowData(sce)%>%
                    data.frame()%>%
                    dplyr::arrange(desc(VAF))%>%
                    dplyr::pull(final_annot),
                  rownames(sce)),]

  print("Adding on protein Data")
  if(protein==TRUE){
    protein_mat <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/layers/read_counts",index=list(NULL,viable_barcodes))
    rownames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ca/id")
    colnames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ra/barcode",index=list(viable_barcodes))
    SingleCellExperiment::altExp(sce, "Protein") <- SingleCellExperiment::SingleCellExperiment(list(Protein=protein_mat))
  }
  
  print("Tabulating cell QC")
  logical_operation <- function(...) Reduce(`&`, ...)
  
  complete_cells<-sce%>%
    {list(.@assays@data$NGT_mask,
          .@assays@data$AF_mask,
          .@assays@data$DP_mask,
          .@assays@data$GQ_mask)}%>%
    logical_operation%>%
    data.frame%>%
    dplyr::select_if(~all(. == TRUE))%>%
    {ifelse(colnames(sce@assays@data$NGT)%in%colnames(.), "Complete", "Other")}
  print(table(complete_cells))
  existing_metadata <- SummarizedExperiment::colData(sce)
  existing_metadata$Required<-complete_cells
  SummarizedExperiment::colData(sce)<-existing_metadata
  
  # Bobby wants this moved to the enumerate clones for some reason instead of just doing it here.
  #print("Computing clones")
  #clone_code<-apply(sce@assays@data$NGT,2,function(x){
  #  paste(ifelse(is.na(x),3,x), sep = "_", collapse = "_")
  #})
 
  sce@metadata$file<-file
  return(sce)
}