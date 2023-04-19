
#' Active learning approach to iteratively filter NGT matrix for missing genotype information.
#'
#' @param NGT_filter required, filtered NGT matrix containing 0,1,2 and NA for cells on rows, and variants in columnns
#' @param variants_of_interest required, universe of potential variants that have been prefiltered
#' @param required_variants optional, vector of variants that cannot be removed
#' @param required_cells optional, vector of cells that cannot be removed
#' @param variant_score_cutoff required, convergence value for fraction of cells that a given variant is successfully genotyped for scale is 0-1
#' @param cell_score_cutoff required,  convergence value for fraction of variants that a given cell is successfully genotyped for scale is 0-1
#' @param greedy_scalar optional, scalar value to speed up computation and remove more cells/variants at each iteration, keep at 0.01 or less
#' @importFrom rlang .data
#' @return NGT_subset, a matrix with reduced dimensions from NGT_filter that possesses well represented variants and cells
#' @export
#'
#' @examples
active_NGT_filter<-function(NGT_filter,
                            variants_of_interest,
                            required_variants,
                            required_cells,
                            variant_score_cutoff=0.9,
                            cell_score_cutoff=0.9,
                            greedy_scalar=0.01){

              kept_variants<-setdiff(colnames(NGT_filter),"Cell")
              kept_cells<-NGT_filter$Cell
              variants_to_remove<-c()
              cells_to_remove<-c()
              variants_to_remove_temp <-data.frame("Variant"="","Score"="")
              cells_to_remove_temp<- data.frame("Cell"="","Score"="")
              NGT_long<-NGT_filter%>%
                            tidyr::pivot_longer(cols=!.data$Cell,names_to="Variant",values_to = "NGT")%>%
                            dplyr::filter(.data$Variant%in%tidyselect::all_of(variants_of_interest))%>%
                            dplyr::mutate(Required_cell=ifelse(.data$Cell%in%tidyselect::all_of(required_cells),TRUE,FALSE))%>%
                            dplyr::mutate(Required_variant=ifelse(.data$Variant%in%tidyselect::all_of(required_variants),TRUE,FALSE))

              while(length(variants_to_remove_temp$Variant)!=0& length(cells_to_remove_temp$Cell)!=0){
                          variants_to_remove_temp<-NGT_long%>%
                                            dplyr::group_by(.data$Variant)%>%
                                            dplyr::summarize(Score=sum(!is.na(.data$NGT))/length(.data$NGT),
                                                                                    required=sum(.data$Required_variant))%>%
                                            dplyr::filter(!.data$Variant%in%tidyselect::all_of(required_variants))%>%
                                            dplyr::filter(.data$Score<=(min(.data$Score)+min(.data$Score)*tidyselect::all_of(greedy_scalar)))%>%
                                            dplyr::filter(.data$Score<tidyselect::all_of(variant_score_cutoff))

                          variants_to_remove<-unique(c(variants_to_remove,variants_to_remove_temp$Variant))

                          cells_to_remove_temp<-NGT_long%>%
                                            dplyr::group_by(.data$Cell)%>%
                                            dplyr::summarize(Score=sum(!is.na(.data$NGT))/length(.data$NGT),
                                                                                required=sum(.data$Required_cell))%>%
                                            dplyr::filter(!.data$Cell%in%tidyselect::all_of(required_cells))%>%
                                            dplyr::filter(.data$Score<=(min(.data$Score)+min(.data$Score)*tidyselect::all_of(greedy_scalar)))%>%
                                            dplyr::filter(.data$Score<tidyselect::all_of(cell_score_cutoff))
                          cells_to_remove<-unique(c(cells_to_remove,cells_to_remove_temp$Cell))

                          NGT_long%<>%
                                          dplyr::filter(!.data$Variant%in%tidyselect::all_of(variants_to_remove))%>%
                                          dplyr::filter(!.data$Cell%in%tidyselect::all_of(cells_to_remove))
                                 }

              NGT_subset<-NGT_long%>%
                            dplyr::filter(!.data$Cell%in%tidyselect::all_of(cells_to_remove))%>%
                            dplyr::filter(!.data$Variant%in%tidyselect::all_of(variants_to_remove))%>%
                            dplyr::select(.data$Cell,.data$Variant,.data$NGT)%>%
                            tidyr::pivot_wider(names_from=.data$Variant,values_from=.data$NGT)

              return(NGT_subset)
}
