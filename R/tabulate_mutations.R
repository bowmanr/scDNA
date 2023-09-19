#' Title
#'
#' @param file 
#' @param GT_cutoff 
#' @param VAF_cutoff 
#'
#' @return
#' @export
#'
#' @examples
tabulate_mutations<-function(file,
                               GT_cutoff=0,
                               VAF_cutoff=0){
  
  file<- "~/CodingCamp/Projects/scDNA_demo/Sample1962.dna+protein.h5"
  
  print(paste("Input file:",file))
  sample_set<-rhdf5::h5read(file=file,name="/assays/dna_variants/metadata/sample_name")[1,]
  
  print(paste("Detected n =",length(sample_set),"sample(s):", paste(sample_set,sep=" ",collapse = " ")))
  
  print("Checking for Barcode Duplicates")
  all_barcodes<- rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
  dedup_barcodes<-all_barcodes[!(duplicated(all_barcodes) | duplicated(all_barcodes, fromLast = TRUE))]
  viable_barcodes<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%dedup_barcodes)
  print(paste(length(all_barcodes)-length(viable_barcodes), "duplicated barcodes detected and removed"))
  
  print("Identifying Variants Past Threshold")
  print(paste("Genotyping Threshold",GT_cutoff))
  print(paste("VAF Threshold",VAF_cutoff, "For Any Sample (if merged)"))
  
  total_variants<-lapply(sample_set,function(sample){
    print(paste("Loading Genotype Data for:",sample))
    sample_of_interest<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/sample_name")%in%sample)
    sample_index<-intersect(sample_of_interest,viable_barcodes)
    
    NGT_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(NULL,sample_index))%>%
      {as(.,"dgCMatrix")}%>%
      `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(sample_index)))  %>% 
      `rownames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")) 
    
    return(data.frame("barcode" = colnames(NGT_data), 
                      "WT"=sparseMatrixStats::colCounts(NGT_data,value = 0),
                      "Het"=sparseMatrixStats::colCounts(NGT_data,value = 1),
                      "Hom"=sparseMatrixStats::colCounts(NGT_data,value = 2),
                      "Missing"=sparseMatrixStats::colCounts(NGT_data,value = 3))%>%sfa
    )
  })
  # The following if else are the changes made to this file to pull in all the variants.
  # I think the else needs to be changed to handle multiple samples correctly?
  # Plus another input needs to be selected that might be missing.
  if(length(sample_set)==1){
    annotated_variants<- scDNA::annotate_variants(file,select_variants = NULL)
    out<-annotated_variants%>% 
      dplyr::full_join(total_variants[[1]],by="id")%>%
      data.frame()%>%
      dplyr::select(id,SYMBOL,AA_change,CONSEQUENCE,WT,Het,Hom,Missing,VAF,genotyping_rate)%>%
      dplyr::arrange(desc(VAF))
    
    return(out)
  } else{ 
    annotated_variants<- scDNA::annotate_variants(file,select_variants = NULL)
    out<-annotated_variants%>% 
      dplyr::full_join(total_variants,by="id")%>%
      data.frame()%>%
      dplyr::select(id,SYMBOL,AA_change,CONSEQUENCE,WT,Het,Hom,Missing,VAF,genotyping_rate)%>%
      dplyr::arrange(desc(VAF))
    
    names(total_variants) <- sample_set
    return(out)
  }
}