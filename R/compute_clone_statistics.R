#' Title
#'
#' @param sce input sce object
#' @param clone_size_cutoff lower 95% confidence interval cutoff for clone inclusion default=10
#'
#' @return
#' @export
#'
#' @examples
compute_clone_statistics<-function(sce,
                                   clone_size_cutoff=10,
                                  skip_ploidy=TRUE){
  
  if(("Group"%in%colnames(sce@metadata$Clones))) {
    stop(message("Clone QC already assessed"))
  }
  
  # Only changes to making them the sce object not the old dataframe format
  print("Computing clone level statistics")
  sce@metadata$Clones<-data.frame(
      sce@assays@data$AF,row.names = rownames(sce))%>%
      tibble::rownames_to_column(var="variants")%>%
      tidyr::pivot_longer(cols=!c(variants),
                            names_to="Cell",
                            values_to="AF")%>%
      mutate(DP=sce@assays@data$DP%>%data.frame()%>%
               tidyr::pivot_longer(cols=tidyselect::everything(),
                                   names_to="Cell",
                                   values_to="DP")%>%
               dplyr::pull(DP))%>%
     mutate(GQ=sce@assays@data$GQ%>%data.frame()%>%
                 tidyr::pivot_longer(cols=tidyselect::everything(),
                                     names_to="Cell",
                                     values_to="GQ")%>%
                 dplyr::pull(GQ))%>%
     mutate(NGT=sce@assays@data$DP%>%data.frame()%>%
                tidyr::pivot_longer(cols=tidyselect::everything(),
                                    names_to="Cell",
                                    values_to="NGT")%>%
                dplyr::pull(NGT)) %>%
    dplyr::inner_join(sce@metadata$NGT%>%
                        dplyr::select(Cell,Group,Clone),by="Cell")%>%
    dplyr::group_by(Group,Clone,variants)%>%
    dplyr::reframe(AF_med=median(AF),
                     DP_med=median(DP),
                     GQ_med=median(GQ))%>%
    dplyr::inner_join(sce@metadata$Clonal_Abundance,by="Clone")
  
  if(sum(sce@metadata$Clones$Group=="Complete")==0){
  print("No cells with complete genotyping, skipping statistcal enumeration")
  } else{
  print("Computing sample level statistics")
  sce@metadata$sample_stats<-   data.frame("Shannon"=sce@metadata$Clones%>%
                                             dplyr::ungroup()%>%
                                             dplyr::filter(Group=="Complete")%>%
                                             dplyr::filter(n_Complete>clone_size_cutoff)%>%
                                             dplyr::distinct(Clone,n_Complete)%>%
                                             dplyr::reframe("Shannon"=vegan::diversity(n_Complete,index="shannon"),.groups=NULL) %>%
                                             dplyr::pull(Shannon),
                                           "Number_of_clones"=sce@metadata$Clones%>%
                                             dplyr::ungroup()%>%
                                             dplyr::filter(Group=="Complete")%>%
                                             dplyr::filter(n_Complete>clone_size_cutoff)%>%
                                             dplyr::distinct(Clone,n_Complete)%>%
                                             dplyr::reframe("Number_of_clones"=nrow(.)) %>%
                                             dplyr::pull(Number_of_clones),
                                           "Number_of_mutations"=length(rownames(sce)),
                                           "Number_of_mutations_in_dominant_clone"=sce@metadata$Clones%>%
                                             dplyr::ungroup()%>%
                                             dplyr::filter(Group=="Complete")%>%
                                             dplyr::filter(grepl("1|2",Clone))%>%
                                             dplyr::distinct(Clone,n_Complete)%>%
                                             dplyr::filter(n_Complete==max(n_Complete))%>%
                                             dplyr::pull(Clone)%>%unique()%>%
                                             strsplit(.,split="_")%>%unlist%>%as.numeric%>%sum,
                                           "Dominant_clone_size"=sce@metadata$Clones%>%
                                             dplyr::ungroup()%>%
                                             dplyr::filter(Group=="Complete")%>%
                                             dplyr::distinct(Clone,n_Complete)%>%
                                             dplyr::mutate(WT_clone=case_when(
                                               !grepl("1|2",Clone)~"WT",
                                               TRUE~"Mutant"))%>%
                                             dplyr::group_by(WT_clone)%>%
                                             dplyr::mutate(Dominant_clone=case_when(
                                               n_Complete==max(n_Complete)~"Dominant",
                                               TRUE~"Other"
                                             ))%>%
                                             dplyr::ungroup()%>%
                                             dplyr::reframe("Dominant_clone_size"=n_Complete[Dominant_clone=="Dominant"&WT_clone!="WT"]/sum(n_Complete))%>%
                                             dplyr::pull("Dominant_clone_size"))
  
      
  }                                             
  if(skip_ploidy){  
    print("Computing Ploidy")   
    sce<-readDNA_CN_H5(sce)
  }
  
  return(sce)
}
