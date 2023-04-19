#' Annotate variants of interest
#'
#' @param file path to h5 file
#' @param annotation_key external data frame containing annotation information
#' @param txdb txdb object for genes of interest
#' @param banned data.frame of genomic coordinates to be filtered out
#' @param NGT filtered NGT matrix from quality_filter_NGT
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom methods as
#' @importFrom rlang .data
#' @return variant annotation matrix
#' @export
#' @examples
annotate_variants<- function(file,annotation_key,txdb,banned,NGT){
    variants<-rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")
    select_variants<-setdiff(colnames(NGT),"Cell")
    SNV_mat<-data.frame(do.call(cbind,
                              rhdf5::h5read(file=file,name="/assays/dna_variants/ca/"))) %>%
                        dplyr::filter(.data$id%in%tidyselect::all_of(select_variants))%>%
                        dplyr::filter(!.data$id%in%tidyselect::all_of(banned))%>%
                        dplyr::mutate(ALT=gsub("\\*","N",.data$ALT))%>%
                        dplyr::mutate(REF=gsub("\\*","N",.data$REF))%>%
                        dplyr::mutate(CHROM=paste0("chr",.data$CHROM))

  #necessary for meaningful GRangess
  SNV_mat$REF<-as(SNV_mat$REF, "DNAStringSet")
  SNV_mat$ALT<-as(SNV_mat$ALT, "DNAStringSet")

  variant_gRange<-GenomicRanges::makeGRangesFromDataFrame(SNV_mat,
                                                           seqnames.field = "CHROM",
                                                           start.field="POS",
                                                           end.field="POS",
                                                           keep.extra.columns=TRUE)

  #necessary for downstream joining of
  variant_gRange$QUERYID<-1:length(variant_gRange)

  #identify and isolate non coding variants
  non_coding_variants <- VariantAnnotation::locateVariants(query=variant_gRange,
                                                           subject=txdb,
                                                           region=VariantAnnotation::AllVariants())%>%
                          data.frame()%>%
                          plyranges::filter(as.character(.data$LOCATION)!="coding")%>%
                          dplyr::inner_join(variant_gRange,by="QUERYID",copy=TRUE)

  #identify and isolate  coding variants
  coding_variants  <-  VariantAnnotation::predictCoding(query=variant_gRange,
                                                        subject=txdb,
                                                        seqSource=Hsapiens,
                                                        varAllele=variant_gRange$ALT)%>%
                        data.frame()

  #Bind it all together into one big table.
  out <- dplyr::bind_rows(non_coding_variants,coding_variants) %>%
                  dplyr::inner_join(annotation_key)%>%
                  dplyr::mutate(AA=ifelse(!is.na(.data$REFAA),
                                   paste0(.data$gene_name,".",.data$REFAA,.data$PROTEINLOC,.data$VARAA),
                                   paste0(.data$gene_name,".intronic")))%>%
                  dplyr::filter(.data$id%in%tidyselect::all_of(colnames(NGT)))

  #append Bulk VAF for reference in future cutoffs and allele selection
  final_mutation_info<-data.frame(
    "Bulk_VAF"=apply(NGT%>%dplyr::select(!.data$Cell),MARGIN=2,function(x){
                                        (sum(x,na.rm=TRUE) / (sum(!is.na(x))*2))*100
                                   }),
   "GT_call_rate"=apply(NGT%>%dplyr::select(!.data$Cell),MARGIN=2,function(x){
                                       100- (sum(is.na(x)) / length(x)*100)
                                    }),
   "id"=colnames(NGT)[colnames(NGT)!="Cell"])%>%
    dplyr::inner_join(out,by="id")

  return(final_mutation_info)
}
