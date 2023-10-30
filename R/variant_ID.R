#' Variant identification and frequency tallies
#'
#' @param file path to the h5 file
#' @param panel name of prebuilt panel/txdb
#' @param GT_cutoff Fraction of cells that are successfully genotyped for initial filtering (default 0.2, meaning 20%)
#' @param VAF_cutoff Fraction of cells that are mutated for initial filtering of variants (default 0.005, meaning 0.05%)
#'
#' @return A dataframe with each variant on a row, and tally of the number of cells that are WT, Het, Hom or missing for a mutation. Calculated VAF and gentoyping freuqency is also provided. If multiple samples are present in the h5 file, a list object will be returned with each sample as an entry in the list
#' @export
#'
#' @examples
variant_ID<-function(file,
                     panel=NULL,
                      GT_cutoff=0,
                      VAF_cutoff=0){

  if(is.null(panel)){
    panel_extract<-rhdf5::h5read(file = file,name = "/assays/dna_read_counts/metadata/panel_name")[1]
    print(paste(panel_extract,"panel detected in H5 file"))
    if(panel_extract%in%c("Myeloid","MSK_RL")){
      print("Using prebuilt TxDB")
      panel<-panel_extract
    } else {
      print("TxDB not detected. To make a panel specific TxDB object use the 'generate_txdb' function")
      print("Defaulting to genome wide UCSC TxDB")
      panel<- "UCSC"
    }
  } else {
    if(panel%in%c("Myeloid","MSK_RL")){
      print("Using prebuilt TxDB")
      panel<-panel
    } else {
      print("TxDB not detected, check spelling. To make a panel specific TxDB object use the 'generate_txdb' function")
      print("Defaulting to genome wide UCSC TxDB")
      panel<- "UCSC"
    }
  }

  if(grepl("loom",file)){
  print(paste("Input file is loom :",file))
    sample_set<- "placeholder"
    print("Identifying Variants Past Threshold")
    print(paste("Genotyping Threshold",GT_cutoff))
    print(paste("VAF Threshold",VAF_cutoff, "For Any Sample (if merged)"))
      
      NGT_data<-rhdf5::h5read(file=file,name="/matrix")%>%
        {as(.,"dgCMatrix")}%>%t()%>%
        `colnames<-`(., paste0("Cell",1:ncol(.))) %>%
        `rownames<-`(., rhdf5::h5read(file=file,name="/row_attrs/id") ) 
      
      total_variants<-list(data.frame("id" = rhdf5::h5read(file=file,name="/row_attrs/id"), 
                        "WT"=sparseMatrixStats::rowCounts(NGT_data,value = 0),
                        "Het"=sparseMatrixStats::rowCounts(NGT_data,value = 1),
                        "Hom"=sparseMatrixStats::rowCounts(NGT_data,value = 2),
                        "Missing"=sparseMatrixStats::rowCounts(NGT_data,value = 3))%>%
               dplyr::mutate(VAF=((Het+Hom*2)/((WT+Het+Hom)*2))*100)%>%
               dplyr::mutate(genotyping_rate=((WT+Het+Hom)/(WT+Het+Hom+Missing) )*100)%>%
               dplyr::filter(genotyping_rate>=GT_cutoff)%>%
               dplyr::filter(VAF>VAF_cutoff)%>%
               dplyr::arrange(desc(VAF)))

  
  } else if(grepl("h5",file)){

  print(paste("Input file is h5 :",file))
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
      
      return(data.frame("id" = rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id"), 
                        "WT"=sparseMatrixStats::rowCounts(NGT_data,value = 0),
                        "Het"=sparseMatrixStats::rowCounts(NGT_data,value = 1),
                        "Hom"=sparseMatrixStats::rowCounts(NGT_data,value = 2),
                        "Missing"=sparseMatrixStats::rowCounts(NGT_data,value = 3))%>%
               dplyr::mutate(VAF=((Het+Hom*2)/((WT+Het+Hom)*2))*100)%>%
               dplyr::mutate(genotyping_rate=((WT+Het+Hom)/(WT+Het+Hom+Missing) )*100)%>%
               dplyr::filter(genotyping_rate>=GT_cutoff)%>%
               dplyr::filter(VAF>VAF_cutoff)%>%
               dplyr::arrange((VAF))
      )
    })
  } 
    # The following if else are the changes made to this file to pull in all the variants.
    # I think the else needs to be changed to handle multiple samples correctly?
    # Plus another input needs to be selected that might be missing.
    if(length(sample_set)==1){ 
        annotated_variants<- annotate_variants(file,panel=panel,select_variants=total_variants$id)
        out<-annotated_variants%>% 
          dplyr::full_join(total_variants[[1]] ,by="id")%>%
          dplyr::filter(!is.na(id))%>%
          dplyr::mutate(final_annot=dplyr::case_when(
            is.na(final_annot)~id,
            TRUE~final_annot))%>%
          data.frame()
        if(nrow(out)==nrow(total_variants[[1]] )){
          print("All variants accounted for")
        } else {
          print(paste("Lost",nrow(total_variants[[1]])-nrow(out),"variants out of",nrow(total_variants[[1]]), "total variants due to poor annotation."))
        }
        return(out)             
    } else if(length(sample_set)==2) {
      total_variants_new<-dplyr::full_join(total_variants[[1]],
                                           total_variants[[2]],
                                          by="id",
                                          suffix=c(paste0("_",sample_set[1]),
                                                    paste0("_",sample_set[2])))
      annotated_variants<- annotate_variants(file,panel=panel,select_variants=total_variants_new$id)
      out<-annotated_variants%>% 
        dplyr::full_join(total_variants_new,by="id")%>%
        dplyr::filter(!is.na(id))%>%
        dplyr::mutate(final_annot=dplyr::case_when(
          is.na(final_annot)~id,
          TRUE~final_annot))%>%
        data.frame()
      if(nrow(out)==nrow(total_variants_new)){
        print("All variants accounted for")
      } else {
        print(paste("Lost",nrow(total_variants_new)-nrow(out),"variants out of",nrow(total_variants[[1]]), "total variants due to poor annotation."))
              }
      return(out)             
  } else if(length(sample_set)>2) {
      print("Not currently functionining for multiple samples")
      print("Returning list of variants without annotation")
      return(total_variants)
    }
 
}

