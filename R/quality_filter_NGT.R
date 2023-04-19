#' Filter NGT matrix based on read depth, allele frequency, and genotype quality
#'
#' @param file path to H5 file
#' @param NGT NGT object from 'read_tapestri_h5_NGT' function
#' @param DP_cut read depth cutoff
#' @param AF_cut allele frequency cutoff
#' @param GQ_cut genotype quality cutoff#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @return filtered NGT matrix containing NAs for cells that do not pass QC
#' @export
#' @examples
quality_filter_NGT<-function(file,NGT,DP_cut=10,AF_cut=20,GQ_cut=20){
  variants<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")
  select_variants<-match(rownames(NGT),variants)

  AF<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/AF",index=list(select_variants,NULL))
  DP<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/DP",index=list(select_variants,NULL))
  GQ<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/GQ",index=list(select_variants,NULL))


  #bind together NGTlim AF, DP, GQ aNGT_mask&ANGT_mask&A
  NGT_filter<-data.frame(setNames(

    # produce long form allele frequency data
    data.frame(AF,
               "variants"=rownames(NGT)),
    c(tidyselect::all_of(colnames(NGT)),"variants")) %>%
      tidyr::pivot_longer(cols=!c(variants),
                   names_to="Cell",
                   values_to="AF"),

    #produce long form allele depth data
    data.frame(DP)%>%
      tidyr::pivot_longer(cols=tidyselect::everything(),
                   names_to="Cell",
                   values_to="DP")%>%dplyr::select(DP),

    #produce long form genotype quality data
    data.frame(GQ)%>%
      tidyr::pivot_longer(cols=tidyselect::everything(),
                   names_to="Cell",
                   values_to="GQ")%>%dplyr::select(GQ),

    #produce long form genotype call data
    data.frame(NGT)%>%
      tidyr::pivot_longer(cols=tidyselect::everything(), names_to="Cell",
                   values_to="NGT")%>%dplyr::select(NGT)) %>%

    #filter DP and GQ
    dplyr::filter(DP>=DP_cut&
             GQ>=GQ_cut)%>%

    #filter AF for each genotype call
    dplyr::mutate(pass=dplyr::case_when(
      NGT==1&(AF>AF_cut)&(AF<(100-AF_cut)) ~ "include",
      NGT==1&((AF<=AF_cut)|(AF>=(100-AF_cut))) ~ "exclude",
      NGT==2&AF>=(100-AF_cut) ~ "include",
      NGT==2&AF<(100-AF_cut) ~ "exclude",
      NGT==0&AF<=AF_cut ~ "include",
      NGT==0&AF>AF_cut ~ "exclude",
      TRUE ~"other"
    ))%>%

    #filter for cells that pass allele frequency cutoff
    dplyr::filter(.data$pass=="include")%>%

    #reform NGT matrix, where missing values are now now
    tidyr::pivot_wider(id_cols=.data$Cell,
                           names_from=.data$variants,
                           values_from=.data$NGT)

return(NGT_filter)
}
