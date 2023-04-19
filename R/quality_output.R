#' Produce long form quality metrics for each cell-variant pair for read depth, allele frequency, and genotype quality
#'
#' @param file path to H5 file
#' @param NGT NGT object from 'read_tapestri_h5_NGT' function
#' @param DP_cut read depth cutoff
#' @param AF_cut allele frequency cutoff
#' @param GQ_cut genotype quality cutoff
#' @param filter default=TRUE, determines whether function should output filtered or unfiltered long form matrix containing per cell-variant pair of NGT, DP, AF, GQ metrics
#' @param input_variants variants of interest to collect QC data
#' @param input_cells cells of interest to subset on
#' @param NGT NGT matrix
#' @importFrom stats setNames
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
#' @importFrom tidyselect everything
#' @importFrom rlang .data
#' @return filtered NGT matrix containing NAs for cells that do not pass QC
#' @export
#'
#' @examples
quality_output<-function(file,
                          filter=TRUE,
                          input_variants,
                          input_cells,
                          NGT,
                          DP_cut=10,
                          AF_cut=20,
                          GQ_cut=20) {
  all_variants<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")
  select_variants<-match(input_variants,all_variants)

  all_cells<-rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
  select_cells<-match(input_cells,all_cells)

  AF<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/AF",index=list(select_variants,select_cells))
  DP<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/DP",index=list(select_variants,select_cells))
  GQ<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/GQ",index=list(select_variants,select_cells))

  ADO<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/ado_rate")
  ADO_select<-which(ADO!="-1")
  NGT_ADO<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(ADO_select,select_cells))

  haplotype<-apply(NGT_ADO,1,median)
  ADO_rate_by_cell<-data.frame(
                          "Cell"=input_cells,
                          "ADO"=apply(NGT_ADO,2,function(x){sum(x!=haplotype)/length(x)}))

  #bind together NGTlim AF, DP, GQ aNGT_mask&ANGT_mask&A
  filtered_long<-data.frame(setNames(

    # produce long form allele frequency data
    data.frame(AF,
               "variants"=input_variants),
    c(tidyselect::all_of(input_cells),"variants")) %>%
      tidyr::pivot_longer(cols=!c(variants),
                              names_to="Cell",
                              values_to="AF"),

    #produce long form allele depth data
    data.frame(DP)%>%
      tidyr::pivot_longer(cols=everything(),
                             # names_to="Cell",
                              values_to="DP")%>%dplyr::select(DP),

    #produce long form genotype quality data
    data.frame(GQ)%>%
      tidyr::pivot_longer(cols=everything(),
                            #  names_to="Cell",
                              values_to="GQ")%>%dplyr::select(GQ),

    #produce long form genotype call data
     NGT%>%dplyr::select(tidyselect::all_of(input_variants))%>%
      tidyr::pivot_longer(cols=everything(),
                            #  names_to="Cell",
                              values_to="NGT")%>%dplyr::select(NGT))%>%
    dplyr::inner_join(ADO_rate_by_cell,by="Cell")


    if(filter==TRUE){
    #filter DP and GQ
     filtered_long%<>%dplyr::filter(DP>=DP_cut&
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
        dplyr::filter(.data$pass=="include")
    }


  return(filtered_long)
}
