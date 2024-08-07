

#' Enumerate clones
#'
#' @param sce  SingleCellExperiment object containing NGT matrix for clone identification.
#' @param replicates number of bootstrapping replicates
#' @importFrom rlang .data
#' @importFrom magrittr %<>%
#' @return updated sce containing a table of clones, modified NGT matrix, and clonal architecture.
#' @export
#'
#' @examples
enumerate_clones<-function(sce,
                           replicates=100
                           ){
  
  print("Computing clones")

  print("Tabulating Cell QC")
  logical_operation <- function(...) Reduce(`&`, ...)
  
  complete_cells <- sce %>% {
    list(.@assays@data$NGT_mask, 
         .@assays@data$AF_mask, 
         .@assays@data$DP_mask, 
         .@assays@data$GQ_mask)
  }%>% 
    logical_operation%>%
    {colSums(.)/nrow(.)}%>%
    {ifelse(.==1, "Complete", "Other")}
  
  print(table(complete_cells))
  existing_metadata <- SummarizedExperiment::colData(sce)
  existing_metadata$Required<-complete_cells
  SummarizedExperiment::colData(sce)<-existing_metadata
  
  reordered_NGT<-sce@assays@data$NGT
  clone_code<-apply(reordered_NGT,2,function(x){
    paste(ifelse(is.na(x),3,x), sep = "_", collapse = "_")
  })
  existing_metadata <- SummarizedExperiment::colData(sce)
  existing_metadata$clone_code<-clone_code
  SummarizedExperiment::colData(sce)<-existing_metadata
  
  #alternatively rename rows with sce@assays@data@metadata$SYMBOL[order(match(sce@assays@data@metadata$annotation_table$id,rownames(sce)))]
  NGT<- data.frame(Cell=colnames(reordered_NGT),
                   Group=SummarizedExperiment::colData(sce)$Required,
                   Clone=SummarizedExperiment::colData(sce)$clone_code,
                   t(reordered_NGT))
  colnames(NGT)<- c("Cell","Group","Clone",rownames(sce))
  
  # NGT to clone is the non-3s while NGT is with the 3's
  NGT_to_clone <- dplyr::filter(NGT,if_all(!c(Cell,Group,Clone),~.< 3))
  
  cols<-c("Clone"=NA,"Complete"=NA,"Other"=NA)
  
  if("Group"%in%colnames(NGT)){
    clonal_abundance <- NGT_to_clone%>%
                              dplyr::group_by(Group)%>%
                              dplyr::count(Clone,name="Count")%>%
                              dplyr::arrange(Count)%>%
                              tidyr::pivot_wider(id_cols=Clone,names_from = Group,values_from = Count)%>%
                              tibble::add_column(!!!cols[!names(cols) %in% names(.)])%>%
                              dplyr::mutate(Complete=ifelse(is.na(Complete),0,Complete))%>%
                              dplyr::mutate(Other=ifelse(is.na(Other),0,Other))%>%
                              dplyr::group_by(Clone)%>%
                              dplyr::mutate(Count=sum(Other,Complete))%>%
                              dplyr::ungroup()%>%
                              dplyr::select(Clone,Clone,Complete,Other,Count)


  } else {
     clonal_abundance <- NGT_to_clone%>%
                          dplyr::count(Clone,name="Count")%>%
                          dplyr::arrange(Count)
  }


# Setup a resampling function to generate multiple clonal abundance tallies
resample_fun<-function(data){
  x <- data[sample(x=1:nrow(data),replace=TRUE),]
  return(as.matrix(x%>%
                      dplyr::count(Clone,name="Count")%>%
                      dplyr::arrange(Count)))
}

print("Computing confidence intervals for all cells")
test<-replicate(n=replicates,resample_fun(NGT_to_clone),simplify = "array")
if(is(test,"list")){
    y <- setNames(lapply(test,data.frame),1:replicates)
  } else if(is(test,"array")){
    y <- setNames(apply(test,3,data.frame),1:replicates)
  }

y <- y%>% purrr::imap(~ purrr::set_names(.x, c("Clone", .y))) %>%
       purrr::reduce(dplyr::full_join, by = "Clone")

y[is.na(y)]<-0

z <- data.frame(t(apply(y%>%dplyr::select(-Clone),1,function(p){
                           stats::quantile(as.numeric(p),probs=c(0.025,0.975))
                    })),
              "Clone"=y$Clone)

clonal_abundance_boot_CI <- setNames(data.frame(
                          dplyr::inner_join(
                                  data.frame(clonal_abundance),
                                   z,
                                    by="Clone")),
                c("Clone","n_Complete","n_Other","Count","LCI","UCI"))


print("Computing confidence intervals for cells with complete genotypes")
NGT_to_clone_complete<-NGT_to_clone%>%dplyr::filter(Group=="Complete")

if(nrow(NGT_to_clone_complete)>0){
test<-replicate(n=replicates,resample_fun(NGT_to_clone_complete),simplify = "array")
if(is(test,"list")){
  y <- setNames(lapply(test,data.frame),1:replicates)
} else if(is(test,"array")){
  y <- setNames(apply(test,3,data.frame),1:replicates)
}

y <- y%>% purrr::imap(~ purrr::set_names(.x, c("Clone", .y))) %>%
  purrr::reduce(dplyr::full_join, by = "Clone")

y[is.na(y)]<-0

z <- data.frame(t(apply(y%>%dplyr::select(-Clone),1,function(p){
  stats::quantile(as.numeric(p),probs=c(0.025,0.975))
})),
"Clone"=y$Clone)

colnames(z)<-c("Complete_LCI","Complete_UCI","Clone")

clonal_abundance_boot_CI <- full_join(clonal_abundance_boot_CI,
                                      z,by="Clone")
} else {
  print("No complete cells. Filling with placeholder 0's")
  
  z<- data.frame("Complete_LCI"=0,
                 "Complete_UCI"=0,
                 "Clone"=clonal_abundance_boot_CI$Clone)
  
  clonal_abundance_boot_CI <- full_join(clonal_abundance_boot_CI,
                                        z,by="Clone")
}

clonal_architecture<-NGT_to_clone%>%
                        dplyr::select(!c(Cell,Group))%>%
                        dplyr::distinct()%>%
                        tidyr::pivot_longer(cols=!Clone,
                                     names_to="final_annot",
                                     values_to="Genotype") %>%
                        dplyr::mutate(Genotype=dplyr::case_when(
                                  Genotype==3~"error",
                                  Genotype==0~"WT",
                                  Genotype==1~"Heterozygous",
                                  Genotype==2~"Homozygous",
                                  TRUE~"error"))#%>%
                        #dplyr::inner_join(sce@assays@data@metadata$annotation_table%>%
                       #                     dplyr::filter(id%in%rownames(sce))%>%
                       #                     dplyr::select(id,AA_change),by="id",relationship="many-to-many")%>%
                      #  dplyr::mutate(AA=factor(AA_change,levels=as.character(unique(AA_change))))
  
  if(any(clonal_architecture$Genotype=="error")){
      "something went wrong"
  }
#Newly added sce metadata needs to be added line by line to avoid overwriting previous metadata
sce@metadata$Clonal_Abundance<-clonal_abundance_boot_CI
sce@metadata$NGT<-NGT_to_clone
sce@metadata$NGT_with_missing<-NGT
sce@metadata$Architecture<-clonal_architecture

  return(sce)
  #return(list("Clones"=clonal_abundance_boot_CI,
  #            "NGT"=NGT_to_clone,
  #            "NGT_with_missing"=NGT,
  #            "Architecture"=clonal_architecture))
}
