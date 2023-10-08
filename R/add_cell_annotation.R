#' Title
#'
#' @param sce 
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
add_cell_annotation<-function(sce,df){
  existing_metadata <- SummarizedExperiment::colData(sce)%>%
    data.frame%>%
    tibble::rownames_to_column(var="Cell")
  if(any(grepl("Cell",colnames(df)))){
    new_metadata  <- df%>%
      dplyr::select(Cell,setdiff(colnames(.),colnames(existing_metadata)))%>%
      dplyr::full_join(existing_metadata,.,by="Cell")
    rownames(new_metadata) <- new_metadata$Cell
    SummarizedExperiment::colData(sce)<-S4Vectors::DataFrame(new_metadata%>%dplyr::select(!Cell))
    return(sce)
  } else{
    print("No 'Cell' column in new meta data")
    return(sce)
    
  }
}
