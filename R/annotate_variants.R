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
<<<<<<< Updated upstream
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

=======
annotate_variants<- function(file,
                             txdb=NULL,
                             select_variants=NULL
                             ){
  
  load("/Users/bowmanrl/Projects/R_packages/scDNA/data/cBioPortal_annotation.rDa",verbose = T)
  
    if(is.null(txdb)){
      print("No TXDB provided, defaulting to complete TxDb.Hsapiens.UCSC.hg19.knownGene")
      custom_txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      genes_found<-genes(custom_txdb)$gene_id
      complete_gene_annotation<-cBioPortal_annotation%>%
        dplyr::filter(entrez_gene_id%in%genes_found)%>%
        dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
        dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
      
    } else if(txdb=="MSK_RL"){
      print("Loading TxDB for Myeloid Clonal Evolution (Levine, MSK)/MSK_RL")
      custom_txdb<-loadDb(system.file('data/MSK_RL_txdb', package = 'scDNA')) # loads in as variable annotation_file?
      #custom_txdb<-loadDb("/Users/bowmanrl/Projects/R_packages/scDNA/data/MSK_RL_txdb")# loads in as variable annotation_file?
      #load(system.file('data/cBioPortal_annotation.rDa', package = 'scDNA')) # loads in as variable annotation_file?
      genes_found<-genes(custom_txdb)$gene_id
      complete_gene_annotation<-cBioPortal_annotation%>%
        dplyr::filter(ensembl_canonical_gene%in%genes_found)%>%
        dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
        dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
      
      } else if(txdb=="UCSC"){
      print("Loading TxDb.Hsapiens.UCSC.hg19.knownGene")
      custom_txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      genes_found<-genes(custom_txdb)$gene_id
      complete_gene_annotation<-cBioPortal_annotation%>%
        dplyr::filter(entrez_gene_id%in%genes_found)%>%
        dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
        dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
      
    } else if(file.exists(txdb)){
      print("Loading custom TxBD")
      custom_txdb<-loadDb(txdb)
      genes_found<-genes(custom_txdb)$gene_id
      
      if(any(grepl("ENSG",genes_found))){
        print("Gene ID appears to be ensemble")
        complete_gene_annotation<-cBioPortal_annotation%>%
          dplyr::filter(entrez_gene_id%in%genes_found)%>%
          dplyr::select(ensembl_canonical_gene,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
          dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
      } else if(!any(grepl("ENSG",genes_found))) {
        print("Trying Entrez Gene ID")
        complete_gene_annotation<-cBioPortal_annotation%>%
          dplyr::filter(entrez_gene_id%in%genes_found)%>%
          dplyr::select(ensembl_canonical_gene,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
          dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
       }
        if(nrow(complete_gene_annotation)>1) {
          print(paste0("Entrez Gene ID extracted. n=",nrow(complete_gene_annotation), " genes detected."))
        } else {
        print("Gene ID could not be detected. Suggest remaking txDB so that ensemble gene names are stored in the gene_id column.")
        print("e.g. ENSG00000122025")
        break
      }
    }

  
  print("Load annotation data")

  print("Extracting Variation Matrix")
  if(grepl("loom",file)){
    SNV_mat_prefilter<-data.frame(do.call(cbind, rhdf5::h5read(file=file,name="/row_attrs/"))) 
  } else if(grepl("h5",file)){
    SNV_mat_prefilter<-data.frame(do.call(cbind, rhdf5::h5read(file=file,name="/assays/dna_variants/ca/"))) 
  }
  
  SNV_mat_prefilter<-  SNV_mat_prefilter%>%
    dplyr::select(id,ALT,REF,CHROM,POS,QUAL,amplicon)%>%
    dplyr::mutate(ALT=gsub("\\*","N",.data$ALT))%>%
    dplyr::mutate(REF=gsub("\\*","N",.data$REF))%>%
    dplyr::mutate(CHROM=paste0("chr",.data$CHROM))%>%
    dplyr::mutate(QUAL=as.numeric(QUAL))%>%
    dplyr::mutate(QUERYID=1:nrow(.))
    
  SNV_mat<-SNV_mat_prefilter%>%  
    dplyr::filter(!grepl("N",REF))
  
  print(paste("Removed",c(nrow(SNV_mat_prefilter)-nrow(SNV_mat)),"variants confounded by upstream deletion"))
  
  if (!is.null(select_variants)) {
    SNV_mat <- SNV_mat %>% dplyr::filter(id %in% select_variants)
  }
  
>>>>>>> Stashed changes
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
