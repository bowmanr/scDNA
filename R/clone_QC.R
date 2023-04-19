#' Title
#'
#' @param final_sample_summary
#' @param file
#'
#' @return
#' @export
#'
#' @examples
clone_QC<-function(final_sample_summary,
                   file){

  if(any(grepl("Group",colnames(final_sample_summary$Clones)))) {
    stop(message("Clone QC already assessed"))
  }

  Cell_QC<-quality_output(file,
                           filter=FALSE,
                           input_variants=setdiff(colnames(final_sample_summary$NGT),c("Cell","Group","Clone")),
                           input_cells=final_sample_summary$NGT$Cell,
                           NGT=final_sample_summary$NGT,
                           DP_cut=10,
                           AF_cut=20,
                           GQ_cut=20)%>%
            dplyr::inner_join(final_sample_summary$NGT%>%
                                      dplyr::select(.data$Cell,.data$Group,.data$Clone),by="Cell")%>%
                                      dplyr::group_by(Group,Clone,variants)%>%
                                      dplyr::summarize(AF_med=median(.data$AF),
                                                       DP_med=median(.data$DP),
                                                       GQ_med=median(.data$GQ),
                                                       ADO_med=median(.data$ADO))%>%
            dplyr::inner_join(final_sample_summary$Clones,by="Clone")
}
