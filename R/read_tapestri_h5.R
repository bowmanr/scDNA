
#' Import Tapestri H5 data
#'
#' @param file path to h5 file
#' @param VAF_cutoff % of cells that are mutated for initial filtering of variants
#' @param GT_cutoff % of cells that are successfully genotyped for initial filtering
#'
#' @return list
#' @export
#'
#' @examples
read_tapestri_h5<-function(file,VAF_cutoff=0.005,GT_cutoff=20){
  NGT<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT")
  NGT[NGT==3]<-NA
  GT_cutoff=GT_cutoff
  VAF_cutoff=VAF_cutoff
  VAF_select<-which(apply(NGT,MARGIN=1,function(x){
    (sum(!is.na(x))/length(x))*100>=GT_cutoff&
      (sum(x,na.rm=TRUE)/(sum(!is.na(x))*2))>=VAF_cutoff
  }))

  AF<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/AF",index=list(VAF_select,NULL))
  DP<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/DP",index=list(VAF_select,NULL))
  GQ<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/GQ",index=list(VAF_select,NULL))
  NGTlim<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(VAF_select,NULL))
  NGTlim[NGTlim==3]<-NA
  variants<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_select))
  cell_barcodes <-rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")
  colnames(NGTlim) <-cell_barcodes
  rownames(NGTlim) <- variants
  return(list(variants,cell_barcodes,NGTlim))
}
