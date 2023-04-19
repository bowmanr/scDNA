#' Function to mimic basic variant and cell cutoffs
#'
#' @param NGT_filter required, filtered NGT matrix containing 0,1,2 and NA for cells on rows, and variants in columnns
#' @param variants_of_interest required, universe of potential variants that have been prefiltered
#' @param required_variants optional, vector of variants that cannot be removed
#' @param required_cells optional, vector of cells that cannot be removed
#' @param variant_score_cutoff required, convergence value for fraction of cells that a given variant is successfully genotyped for scale is 0-1
#' @param cell_score_cutoff required,  convergence value for fraction of variants that a given cell is successfully genotyped for scale is 0-1
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
#' @return NGT matrix
#' @export
#'
#' @examples
TI_mimic_filter<- function(NGT_filter,
                           variants_of_interest,
                           required_variants,
                           required_cells,
                           variant_score_cutoff=50,
                           cell_score_cutoff=50) {

                NGT_long<-NGT_filter%>%
                      tidyr::pivot_longer(cols=!.data$Cell,names_to="Variant",values_to = "NGT")%>%
                      dplyr::filter(.data$Variant%in%tidyselect::all_of(variants_of_interest))%>%
                      dplyr::mutate(Required_cell=ifelse(.data$Cell%in%tidyselect::all_of(required_cells),TRUE,FALSE))%>%
                      dplyr::mutate(Required_variant=ifelse(.data$Variant%in%tidyselect::all_of(required_variants),TRUE,FALSE))

                # here we filter for variants that contain genotype information for at least 50% of cells
                final_variants<-NGT_long%>%
                      dplyr::group_by(.data$variants)%>%
                      dplyr::summarize(gt.mv=(length(.data$NGT)/length(tidyselect::all_of(colnames(NGT_filter))))*100)%>%
                      dplyr::filter(.data$gt.mv>=tidyselect::all_of(variant_score_cutoff))%>%
                      dplyr::pull(.data$variants)

                # here we filter for cells that now contain genotype information for at least 50% of our curated set of variants.
                final_cells<-NGT_long%>%
                      dplyr::filter(.data$variants%in%tidyselect::all_of(final_variants))%>%
                      dplyr::group_by(.data$Cell)%>%
                      dplyr::summarize(gt.mc=(length(.data$NGT)/length(tidyselect::all_of(final_variants)))*100)%>%
                      dplyr:: filter(.data$gt.mc>=tidyselect::all_of(cell_score_cutoff))%>%
                      dplyr::pull(.data$Cell)

# lastly we reconstruct a new NGT matrix of cell-genotype pairs that passed the above filters.
                  final_NGT<- NGT_filter%>%
                      dplyr::filter(.data$Cell%in%tidyselect::all_of(final_cells)&
                                    .data$variants%in%tidyselect::all_of(final_variants))%>%
                      tidyr::pivot_wider(id_cols=.data$Cell,
                                  names_from=.data$variants,
                                  values_from=.data$NGT)
}
