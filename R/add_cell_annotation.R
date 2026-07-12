#' Add Cell Annotation to Metadata
#'
#' @description
#' A quick helper function that adds additional metadata to the existing sce. This is primarily used for demultiplexing and obtaining the cluster ID.
#'
#' @param sce the existing single cell experiment
#' @param df The new metadata dataframe wanted to be added. The df needs to have a column named "Cell" which is the cell names (colnames of the sce)
#'
#' @return it returns an sce object with the new updated metadata integrated.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' sce<-add_cell_annotation(sce,umap_df_final)
#' }
#'
add_cell_annotation<-function(sce,df){

  if(any(grepl("Cell",colnames(df)))){
    existing_metadata <- SummarizedExperiment::colData(sce)%>%
      data.frame%>%
      tibble::rownames_to_column(var="Cell")
    new_cols <- setdiff(colnames(df), colnames(existing_metadata))
    df_to_join <- df[, c("Cell", new_cols), drop = FALSE]
    new_metadata <- dplyr::full_join(existing_metadata, df_to_join, by = "Cell")

    rownames(new_metadata) <- new_metadata$Cell
    SummarizedExperiment::colData(sce)<-S4Vectors::DataFrame(new_metadata%>%dplyr::select(!Cell))
  }
  else{
    warning("No 'Cell' column in new meta data")
  }
  return(sce)

}
