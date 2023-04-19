
#' Select clones of interest on the basis QC metrics
#'
#' @param final_sample_summary
#' @param ADO_cut maximum median allele dropout
#' @param GQ_cut mimimum median gene quality score
#' @param DP_cut minimum median sequencing depth
#' @param select_exact default=FALSE, can supply a vector of clones of interest to bypass QC filtering
#' @return final_sample_summary
#' @export
#'
#' @examples
select_clones<- function(final_sample_summary=final_sample_summary,
                         ADO_cut=0.10,
                         GQ_cut=30,
                         DP_cut=10,
                         select_exact=FALSE){

  if(select_exact==FALSE){
  select_clones<- final_sample_summary$Clones%>%
                            dplyr::filter(ADO_med<ADO_cut)%>%
                            dplyr::filter(GQ_med>GQ_cut)%>%
                            dplyr::filter(DP_med>DP_cut)%>%
                            dplyr::group_by(Clone)%>%
                            dplyr::filter(Group%in%c("Complete"))%>%
                            dplyr::pull(Clone)%>%unique()
  } else{
    select_clones<-select_clones
  }
  final_sample_summary$NGT%<>%dplyr::filter(Clone%in%select_clones)
  final_sample_summary$Clones%<>%dplyr::filter(Clone%in%select_clones)
  final_sample_summary$Architecture%<>%dplyr::filter(Clone%in%select_clones)
  return(final_sample_summary)
}
