
#' Import protein data
#'
#' @param file path to h5 file
#'
#' @return
#' @export
#'
#' @examples
read_tapestri_h5_protein<-function(file=file){
          protein_mat<-rhdf5::h5read(file=file,name="/assays/protein_read_counts/layers/read_counts")
          rownames(protein_mat) <- rhdf5::h5read(file=file,name="/assays/protein_read_counts/ca/id")
          colnames(protein_mat)<-  rhdf5::h5read(file=file,name="/assays/protein_read_counts/ra/barcode")

          protein_mat_final<-data.frame("Cell"=colnames(protein_mat),
                                               t(protein_mat))%>%
                             dplyr::mutate(Cell=gsub("-1","",Cell))%>%
                             dplyr::filter(Cell%in%tidyselect::all_of(final_NGT$Cell))
          return(protein_mat_final)
}

