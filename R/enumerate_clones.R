

#' Enumerate clones
#'
#' @param NGT  NGT matrix for clone identification. Should lack NA's
#' @param replicates number of bootstrapping replicates
#' @param variant_metadata metadata
#' @importFrom rlang .data
#' @importFrom magrittr %<>%
#' @return list containing a table of clones, modified NGT matrix, and clonal architecture.
#' @export
#'
#' @examples
enumerate_clones<-function(NGT,
                           variant_metadata,
                           replicates=100
                           ){

  bulk_VAF_order   <- NGT%>%
                           dplyr::select(!c(.data$Cell,.data$Group))%>%
                           colSums%>%
                           sort(decreasing=TRUE)%>%
                           names

  NGT_to_clone     <- NGT[,c("Cell","Group",bulk_VAF_order)]%>%
                                  tidyr::unite("Clone",
                                               tidyselect::all_of(bulk_VAF_order),
                                               sep="_", remove = FALSE)

  if("Group"%in%colnames(final_NGT)){
    clonal_abundance <- NGT_to_clone%>%
                              dplyr::group_by(Group)%>%
                              dplyr::count(.data$Clone,name="Count")%>%
                              dplyr::arrange(.data$Count)%>%
                              tidyr::pivot_wider(id_cols =.data$Clone,names_from = .data$Group,values_from = .data$Count )%>%
                              dplyr::mutate(Complete=ifelse(is.na(.data$Complete),0,.data$Complete))%>%
                              dplyr::mutate(Other=ifelse(is.na(.data$Other),0,.data$Other))%>%
                              dplyr::group_by(.data$Clone)%>%
                              dplyr::mutate(Count=sum(.data$Other,.data$Complete))%>%
                              dplyr::ungroup()

  } else {
     clonal_abundance <- NGT_to_clone%>%
                          dplyr::count(.data$Clone,name="Count")%>%
                          dplyr::arrange(.data$Count)
  }

# Setup a resampling function to generate multiple clonal abundance tallies
resample_fun<-function(data){
  x <- data[sample(x=1:nrow(data),replace=TRUE),]
  return(as.matrix(x%>%
                      dplyr::count(.data$Clone,name="Count")%>%
                      dplyr::arrange(.data$Count)))
}

test<-replicate(n=replicates,resample_fun(NGT_to_clone),simplify = "array")
if(class(test)=="list"){
    y <- setNames(lapply(test,data.frame),1:replicates)
  } else if(class(test)=="array"){
    y <- setNames(apply(test,3,data.frame),1:replicates)
  }

y %<>% purrr::imap(~ set_names(.x, c("Clone", .y))) %>%
       purrr::reduce(dplyr::full_join, by = "Clone")

y[is.na(y)]<-0

z <- data.frame(t(apply(y%>%dplyr::select(-.data$Clone),1,function(p){
                           stats::quantile(as.numeric(p),probs=c(0.025,0.975))
                    })),
              "Clone"=y$Clone)

clonal_abundance_boot_CI <- setNames(data.frame(
                          dplyr::inner_join(
                                  data.frame(clonal_abundance),
                                   z,
                                    by="Clone")),
                c("Clone","n_Complete","n_Other","Count","LCI","UCI"))



clonal_architecture<-NGT_to_clone%>%
                        dplyr::select(!c(.data$Cell,.data$Group))%>%
                        dplyr::distinct()%>%
                        tidyr::pivot_longer(cols=!.data$Clone,
                                     names_to="id",
                                     values_to="Genotype") %>%
                        dplyr::mutate(Genotype=dplyr::case_when(
                                  Genotype==3~"error",
                                  Genotype==0~"WT",
                                  Genotype==1~"Heterozygous",
                                  Genotype==2~"Homozygous",
                                  TRUE~"error"))%>%
                        dplyr::inner_join(variant_metadata%>%dplyr::select(.data$id,.data$AA),by="id")%>%
                        dplyr::mutate(AA=factor(AA,levels=as.character(final_variant_info$AA)))


  if(any(clonal_architecture$Genotype=="error")){
      "something went wrong"
  }

  return(list("Clones"=clonal_abundance_boot_CI,
              "NGT"=NGT_to_clone,
              "Architecture"=clonal_architecture))
}
