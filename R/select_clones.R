
#' Select clones of interest on the basis QC metrics
#'
#' @param sce a SingleCellExperiment object containing the Clones,Architecture, and NGT matrix
#' @param ADO_cut maximum median allele dropout
#' @param GQ_cut mimimum median gene quality score
#' @param DP_cut minimum median sequencing depth
#' @param select_exact default=FALSE, can supply a vector of clones of interest to bypass QC filtering
#' @return a subset of the sce object
#' @export
#'
#' @examples
select_clones<- function(sce,
                         ADO_cut=0.10,
                         GQ_cut=30,
                         DP_cut=10,
                         select_exact=FALSE){
# Only change made was converting the old dataframe to the sce.
  if(select_exact==FALSE){
  select_clones<- sce@metadata$Clones%>%
                            #dplyr::filter(ADO_med<ADO_cut)%>%
                            dplyr::filter(GQ_med>GQ_cut)%>%
                            dplyr::filter(DP_med>DP_cut)%>%
                            dplyr::group_by(Clone)%>%
                            dplyr::filter(Group%in%c("Complete"))%>%
                            dplyr::pull(Clone)%>%unique()
  } else{
    select_clones<-select_exact
  }
  sce@metadata$NGT<-sce@metadata$NGT%>%dplyr::filter(Clone%in%select_clones)
  sce@metadata$Clones<-sce@metadata$Clones%>%dplyr::filter(Clone%in%select_clones)
  sce@metadata$Architecture<-sce@metadata$Architecture%>%dplyr::filter(Clone%in%select_clones)
  return(sce)
}
