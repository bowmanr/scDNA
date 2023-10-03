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
                                   clone_size_cutoff=10){
  
  if(("Group"%in%colnames(sce@metadata$Clones))) {
    stop(message("Clone QC already assessed"))
  }
  
  # Only changes to making them the sce object not the old dataframe format
  print("Computing clone level statistics")
  sce@metadata$Clones<-(quality_output(file=as.character(sce@metadata$file),
                                                              filter=FALSE,
                                                              input_variants=SummarizedExperiment::rowData(sce)$id,
                                                              input_cells= sce@metadata$NGT$Cell,
                                                              NGT=  sce@metadata$NGT,
                                                              DP_cut=10,
                                                              AF_cut=20,
                                                              GQ_cut=20)%>%
    dplyr::inner_join(sce@metadata$NGT%>%
                        dplyr::select(.data$Cell,.data$Group,.data$Clone),by="Cell")%>%
    dplyr::group_by(Group,Clone,variants)%>%
    dplyr::reframe(AF_med=median(AF),
                     DP_med=median(DP),
                     GQ_med=median(GQ),
                     ADO_med=median(ADO))%>%
    dplyr::inner_join(sce@metadata$Clones,by="Clone"))
  
  print("Computing sample level statistics")
  sce@metadata$sample_stats<- data.frame("Shannon"=sce@metadata$Clones%>%
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
                                                   dplyr::mutate(Dominant_clone=case_when(
                                                     WT_clone=="Mutant"&n_Complete==max(n_Complete)~"Dominant",
                                                     TRUE~"Other"))%>%
                                                   dplyr::reframe("Dominant_clone_size"=n_Complete[Dominant_clone=="Dominant"]/sum(n_Complete))%>%
                                                   dplyr::pull("Dominant_clone_size")
          )
                                                    
    print("Computing Ploidy")   
    sce<-readDNA_CN_H5(sce)

  
  return(sce)
}
