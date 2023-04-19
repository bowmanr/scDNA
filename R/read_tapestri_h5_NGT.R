
#' Import Tapestri H5 data and extract genotype matrix
#'
#' @param file path to h5 file
#' @param VAF_cutoff % of cells that are mutated for initial filtering of variants
#' @param GT_cutoff % of cells that are successfully genotyped for initial filtering
#'
#' @return list
#' @export
#'
#' @examples
read_tapestri_h5_NGT<-function(file,VAF_cutoff=0.005,GT_cutoff=20){
  NGT<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT")
  NGT[NGT==3]<-NA
  VAF_select<-which(apply(NGT,MARGIN=1,function(x){
    (sum(!is.na(x))/length(x))*100>=GT_cutoff&
      (sum(x,na.rm=TRUE)/(sum(!is.na(x))*2))>=VAF_cutoff
  }))
  NGTlim<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(VAF_select,NULL))
  NGTlim[NGTlim==3]<-NA
  rownames(NGTlim)<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_select))
  colnames(NGTlim) <-rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
  return(NGTlim)
}
