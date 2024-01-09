#' This function generates the Copy Number by determining the ploidy of each mutation
#'
#' @param sce the single cell experiment object
#' @param reference_cells set cells that will be used to determine allele dropout and false postive rates 
#'
#' @return
#' @export
#'
#' @examples
readDNA_CN_H5<-function(sce,reference_cells=NULL){
  
  # mark reference cells just take them all for now and set it as a variable above.
  # in the future this could be T cells or whatever
  CNV_sce<-SingleCellExperiment::altExp(sce,"CNV")
  if(is.null(reference_cells)){
    reference_cells = colnames(CNV_sce@assays@data$CNV)
  }
  all_amplicon_data<-CNV_sce@assays@data$CNV%>%
    dplyr::mutate(amplicon = rownames(.))%>%
    tidyr::pivot_longer(cols=!amplicon,values_to = "depth",names_to = "barcode")%>%
    dplyr::inner_join(rhdf5::h5read(file=as.character(sce@metadata$file),name="/assays/dna_read_counts/ca",index=list(NULL,reference_cells))%>%
                      base::do.call(cbind,.)%>%data.frame%>%
                      dplyr::rename(amplicon=id),
                      by="amplicon")%>%
    dplyr::mutate(cell_subset=ifelse(barcode%in%reference_cells,"Ref","Other"))%>%
    dplyr::group_by(barcode)%>%
    dplyr::mutate(fraction=depth/sum(depth))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(amplicon)%>%
    dplyr::mutate(norm_frac=(fraction/median(fraction)))%>%
    dplyr::mutate(ploidy=log2(norm_frac/median(norm_frac[cell_subset=="Ref"])))%>%
    dplyr::arrange(CHROM,start_pos)%>%
    #dplyr::ungroup()%>%
    dplyr::mutate(amplicon=factor(amplicon,levels=c(unique(amplicon))))

  
  # associated amplicons with variant information
  mutant_subset_amplicon_data<-dplyr::inner_join(dplyr::inner_join(SummarizedExperiment::rowData(sce)%>%
                                                       data.frame%>%
                                                       dplyr::select(id,amplicon,final_annot),
                                                     sce@metadata$NGT%>%
                                                       tidyr::pivot_longer(cols = !c(Cell,Group,Clone),names_to = "final_annot",values_to = "NGT"),
                                               by="final_annot")%>%
                                              dplyr::rename(barcode=Cell),
                                          all_amplicon_data,
                                            by=c("barcode","amplicon"))%>%
                                dplyr::mutate(final_annot=factor(final_annot,
                                                           levels=SummarizedExperiment::rowData(sce)%>%
                                                                   data.frame%>%
                                                                   dplyr::pull(final_annot)))
  CNV_sce@metadata$full_ploidy<-(all_amplicon_data)
  CNV_sce@metadata$ploidy<-(mutant_subset_amplicon_data)
  sce@metadata$ploidy<-(mutant_subset_amplicon_data)

# Fraction of points outside of confidence interval lower bound over total number of calls is the false positive rate.
  false_positive_cost<-mutant_subset_amplicon_data%>%
    data.frame%>%
    dplyr::group_by(final_annot,NGT)%>%
    dplyr::select(final_annot,NGT,ploidy)%>%
    dplyr::filter(NGT==1)%>%
    dplyr::mutate(mean=mean(ploidy))%>%
    dplyr::mutate(n=n())%>%
    dplyr::mutate(std = sd(ploidy))%>%
    dplyr::mutate(lower_lim95 = mean -qt(0.9998,df=n-1)*std/sqrt(n))%>% # change this rate to give people multiple CI's
    dplyr::ungroup()%>%
    data.frame%>%
    dplyr::select(final_annot,lower_lim95)%>%
    distinct

 
 fp_cost0<-mutant_subset_amplicon_data%>%
   data.frame%>%
   dplyr::group_by(final_annot,NGT)%>%
   dplyr::select(final_annot,NGT,ploidy)%>%
   dplyr::filter(NGT==0)%>%
   dplyr::inner_join(false_positive_cost,by="final_annot")%>%
   mutate(count=ifelse(ploidy<lower_lim95,1,0))%>%
   mutate(summed=sum(count))%>%
   mutate(false_positiveWT = sum(count)/n())%>%ungroup()%>%dplyr::select(final_annot,false_positiveWT)%>%distinct
 
  fp_cost2<-mutant_subset_amplicon_data%>%
    data.frame%>%
    dplyr::group_by(final_annot,NGT)%>%
    dplyr::select(final_annot,NGT,ploidy)%>%
    dplyr::filter(NGT==2)%>%
    dplyr::inner_join(false_positive_cost,by="final_annot")%>%
    dplyr::mutate(count=ifelse(ploidy<lower_lim95,1,0))%>%
    dplyr::mutate(false_positiveHom = sum(count)/n())%>%ungroup()%>%dplyr::select(final_annot,false_positiveHom)%>%distinct

  sce@metadata$FalsePositive <- dplyr::full_join(fp_cost0,fp_cost2,by="final_annot")
  sce@metadata$FalsePositive[is.na(sce@metadata$FalsePositive)]<- 0
  SingleCellExperiment::altExp(sce, "CNV") <- CNV_sce

  return(sce)
}
