
#' Extract genotypes for missing data in filtered NGT matrices
#'
#' @param NGT_to_fill Matrix with missing values that did not pass QC
#' @param perfect_NGT perfect matrix for indexing groups oof cells
#' @param file path to H5 file containing complete NGT matrix.
#' @importFrom rlang .data
#' @importFrom magrittr %<>%
#' @return imperfect NGT matrix with cells and variants of interest
#' @export
#'
#' @examples
generate_lowQC_matrix<-function(NGT_to_fill,
                                complete_NGT,
                                file){

  variants<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")
  cells <-rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
  variant_select <- match(setdiff(colnames(NGT_to_fill),"Cell"),variants)
  cell_select <-  match(NGT_to_fill%>%dplyr::pull(.data$Cell),cells)
  full_NGT<-setNames(data.frame("Cell"=NGT_to_fill%>%dplyr::pull(.data$Cell),
                    t(rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(variant_select,cell_select)))),
                    colnames(NGT_to_fill))
  full_NGT[full_NGT==3]<-NA
  full_NGT%<>%filter(across(.cols = !.data$Cell,.fns = ~ !is.na(.x)))


  full_NGT%<>%dplyr::mutate("Group"=dplyr::case_when(
                                  .data$Cell%in%tidyselect::all_of(complete_NGT$Cell)~"Complete",
                                  TRUE~"Other"))
 return(full_NGT)
}
