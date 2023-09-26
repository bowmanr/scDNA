#' Title
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
readDNA_CN_H5<-function(sce,reference_cells){
  
  
  
  all_amplicon_data<-
    
    #start from sce and pivot_longer the whole mat
    # originally I had amplicon column, it is not stored in row names
    
    #store this as rowData rhdf5::h5read(file=file,name="/assays/dna_read_counts/ca/",index=list(NULL,viable_barcodes))

    tidyr::pivot_longer(cols=!amplicon,values_to = "depth",names_to = "barcode")%>%
    #inner_join rowData with  depth
    dplyr::inner_join(rhdf5::h5read(file=file,name="/assays/dna_read_counts/ca/",index=list(NULL,viable_barcodes))%>%
                      do.call(cbind,.)%>%data.frame%>%
                      dplyr::rename(amplicon=id),
                      by="amplicon")%>%
    
    # mark reference cells just take them all for now and set it as a variable above.
    # in the future this could be T cells or whatever
    mutate(cell_subset=ifelse(barcode%in%reference_cells,"Ref","Other"))%>%
    group_by(barcode)%>%
    # calculate the fraction
    mutate(fraction=depth/sum(depth))%>%
    ungroup()%>%
    
    group_by(amplicon)%>%
    mutate(norm_frac=(fraction/median(fraction)))%>%
    
    mutate(ploidy=log2(norm_frac/median(norm_frac[cell_subset=="Ref"])))%>%
    arrange(CHROM,start_pos)%>%
    mutate(amplicon=factor(amplicon,levels=c(unique(amplicon))))
  
  #associated amplicons with variant information
  mutant_subset_amplicon_data<-inner_join(inner_join(rowData(sce)%>%
                                                       data.frame%>%
                                                       dplyr::select(id,amplicon,final_annot),
                                                     sce@metadata$NGT%>%
                                                       tidyr::pivot_longer(cols = !c(Cell,Group,Clone),names_to = "final_annot",values_to = "NGT"),
                                               by="final_annot")%>%
                                              dplyr::rename(barcode=Cell),
                                          all_amplicon_data,
                                            by=c("barcode","amplicon"))%>%
                                 mutate(final_annot=factor(final_annot,
                                                           levels=rowData(sce)%>%
                                                                   data.frame%>%
                                                                   dplyr::pull(final_annot)))

  ggplot(mutant_subset_amplicon_data,aes(x=as.factor(Clone),y=ploidy,fill=as.factor(NGT)))+
    geom_boxplot()+
    facet_grid(Group~final_annot)
  
  
}
