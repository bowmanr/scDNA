#' Import Tapestri H5 data and extract genotype matrix
#'
#' @param file path to the h5 file
#' @param GT_cutoff Fraction of cells that are successfully genotyped for initial filtering (default 0.2, meaning 20%)
#' @param VAF_cutoff Fraction of cells that are mutated for initial filtering of variants (default 0.005, meaning 0.05%)
#' @param DP_cutoff minimum number of reads necessary for a reliable genotype call in a single cell (default: 10)
#' @param GQ_cutoff minimum genotype quality necessary for a reliable genotype call in a single cell (default: 30)
#' @param AF_cutoff Deviation from 0, 50, or 100% for a reliable call of WT, Het or Hom respectively (default: 25)
#' @param variant_set character vector of variants to be included in the format output by mission bio.
#' @param demultiplex_cells is the cell names which we want to pull, often left NULL for majority of analysis
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
                    demultiplex_cells=NULL,
                    protein=TRUE){
  
  
  if(is.null(variant_set)){
    print("No variants provided, run 'variant_ID' first")
  }
  
  if(!is.null(variant_set)){
    # variant_set has all info, so we need to pull out the id for the h5 file
    VAF_cut_variants<- variant_set%>%pull(id)
    
    print(paste("Loading n=",length(VAF_cut_variants),"variants")) 

    print(paste("Input file:",file))
    sample_set<-rhdf5::h5read(file=file,name="/assays/dna_variants/metadata/sample_name")[1,]
    
    print(paste("Detected n =",length(sample_set),"sample(s):", paste(sample_set,sep=" ",collapse = " ")))
    
    print("Checking for Barcode Duplicates")
    all_barcodes<- rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
    dedup_barcodes<-all_barcodes[!(duplicated(all_barcodes) | duplicated(all_barcodes, fromLast = TRUE))]
    viable_barcodes<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%dedup_barcodes)
    print(paste(length(all_barcodes)-length(viable_barcodes), "duplicated barcodes detected and removed"))

    if(is.null(demultiplex_cells)){
      sample_of_interest<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/sample_name")%in%sample_set)
    }else{
      sample_of_interest<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%demultiplex_cells)
    }
    viable_barcodes<-intersect(sample_of_interest,viable_barcodes)
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
  rownames(sce)<-SummarizedExperiment::rowData(sce)$final_annot

  print("Tabulating Variant QC")
  logical_operation <- function(...) Reduce(`&`, ...)
  variant_genotype_QC <- sce %>% {
                            list(.@assays@data$NGT_mask, 
                                 .@assays@data$AF_mask, 
                                 .@assays@data$DP_mask, 
                                 .@assays@data$GQ_mask)
                                }%>% 
                            logical_operation %>% 
                            {rowSums(.)/ncol(.)*100}
  existing_rowData <- SummarizedExperiment::rowData(sce)
  existing_rowData$variant_QC<-variant_genotype_QC
  SummarizedExperiment::rowData(sce)<-existing_rowData
  
  print("Reordering sce rows (variants) based on bulk VAF")
  # this is newly added so we can get correct order for annotation later.
  sce<- sce[match(SummarizedExperiment::rowData(sce)%>%
                    data.frame()%>%
                    dplyr::arrange(desc(pick(starts_with("VAF"))))%>%
                    dplyr::pull(final_annot),
                  rownames(sce)),]

 
  #### protein starts here
  if(protein==TRUE){
          skip <- TRUE
          skip <- tryCatch( rhdf5::h5read(file = file, 
                                                    name = "/assays/protein_read_counts/layers/read_counts",
                                                    index=list(NULL,viable_barcodes))%>%
                              nrow()%>%{.>0}, 
                      error = function(e) { 
                        print(paste("'Protein' dataset not found"))
                        return(FALSE)
                        })
          if(skip){
            print("Adding Protein data")
            protein_mat <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/layers/read_counts",index=list(NULL,viable_barcodes))
            rownames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ca/id")
            colnames(protein_mat) <- gsub("-","\\.",colnames(protein_mat))
            colnames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ra/barcode",index=list(viable_barcodes))
            SingleCellExperiment::altExp(sce, "Protein") <- SingleCellExperiment::SingleCellExperiment(list(Protein=protein_mat))
          } 
  }
  
  print("Adding Copy Number data")
    amplicon_data<-rhdf5::h5read(file=file,name="/assays/dna_read_counts/layers/read_counts",index=list(NULL,viable_barcodes))%>% data.frame()
    colnames(amplicon_data) <- rhdf5::h5read(file=file,name="/assays/dna_read_counts/ra/barcode",index=list(viable_barcodes))
    colnames(amplicon_data) <- gsub("-","\\.",colnames(amplicon_data))
    rownames(amplicon_data) <- rhdf5::h5read(file=file,name="/assays/dna_read_counts/ca/id",index=list(NULL))
    SingleCellExperiment::altExp(sce, "CNV") <- SingleCellExperiment::SingleCellExperiment(list(CNV=amplicon_data))
  
  # Moved this to enumerate_clones  
  # print("Tabulating cell QC")
  #   logical_operation <- function(...) Reduce(`&`, ...)
  #   complete_cells<-sce%>%
  #     {list(.@assays@data$NGT_mask,
  #           .@assays@data$AF_mask,
  #           .@assays@data$DP_mask,
  #           .@assays@data$GQ_mask)}%>%
  #     logical_operation%>%
  #     data.frame%>%
  #     dplyr::select_if(~all(. == TRUE))%>%
  #     {ifelse(colnames(sce@assays@data$NGT)%in%colnames(.), "Complete", "Other")}
  #   print(table(complete_cells))
  #   existing_metadata <- SummarizedExperiment::colData(sce)
  #   existing_metadata$Required<-complete_cells
  #   SummarizedExperiment::colData(sce)<-existing_metadata
  #  
  sce@metadata$file<-file
  #sce <-readDNA_CN_H5(sce,reference_cells = NULL) # need to move this to after enumerate clones so we can use NGT
  
  return(sce)
}
