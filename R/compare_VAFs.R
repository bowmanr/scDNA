
#' Title
#'
#' @param sce 
#' @param by 
#' @param variants_to_test 
#'
#' @return
#' @export
#'
#' @examples
compare_VAFs<-function(sce,by=NULL,variants_to_test){
  
  if(is.null(by)){
    print("Provide a column in metadata to group by")
    break
  }
  if(is.null(variants_to_test)){
    print("No variant matrix provided")
    break
  }
  
  file<- sce@metadata$file
  variantID_to_test<-variants_to_test%>%pull(id)%>%unique
  variant_index<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")%in%variantID_to_test)
  NGT_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(variant_index,NULL)) %>%
    `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode"))  %>% 
    `rownames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(variant_index))) 
  
  metadata<-SummarizedExperiment::colData(sce)%>%data.frame%>%tibble::rownames_to_column(var="Cell")
  
  if(!by%in%colnames(metadata)){
    print("Column not found in metadata")
    break
  }
  test_col <-colnames(metadata)[colnames(metadata)==by]
  NGT_split<-NGT_data%>%data.frame()%>%
    tibble::rownames_to_column(var="id")%>%
    tidyr::pivot_longer(cols=-id,
                        values_to="NGT",
                        names_to="Cell")%>%
    dplyr::full_join(metadata,by="Cell")%>%
    dplyr::rename(test=dplyr::all_of(test_col))%>%
    split(.,f=.$test) 
  
  NGT_VAFs<-lapply(names(NGT_split),function(group){
    
    NGT_split[[group]]%>%group_by(id)%>%
      dplyr::reframe(WT=sum(NGT==0),
                     Het=sum(NGT==1),
                     Hom=sum(NGT==2),
                     Missing=sum(NGT==3))%>%
      dplyr::mutate(VAF=((Het+Hom*2)/((WT+Het+Hom)*2))*100)%>%
      dplyr::mutate(genotyping_rate=((WT+Het+Hom)/(WT+Het+Hom+Missing) )*100)%>%
      dplyr::arrange(desc(VAF))%>%
      dplyr::full_join(.,variants_to_test%>%dplyr::select(id,final_annot),by=c("id"="id"))%>%
      dplyr::mutate(Group=group)
  })
  return(Reduce(bind_rows,NGT_VAFs))
}