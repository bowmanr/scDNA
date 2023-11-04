



c#' Annotate variants of interest
#' This function takes in h5 files to extract DNA information through given TXDB files (primarily hg19)
#' The gene names are mapped for their specific nucleotide positions. The starting and ending positions
#' for the chromosome are listed. The consequence of the mutation such as synonymous, nonsynonymous,
#' as well as if this is a coding region, on the exon boundry or intronic. Short amino acid changes are also labeled.
#' The variance matrix also extracts the amplicon that the variant is found on.
#' @param file path to h5 file to pull out relevant DNA information
#' @param select_variants variants of interest, default is all variants
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom methods as
#' @importFrom rlang .data
#' @importFrom dplyr n
#' @return variant annotation matrix
#' @export
#' @examples
annotate_variants<- function(file,
                             panel=NULL,
                             select_variants=NULL){
  
  load(system.file(paste0('data/cBioPortal_annotation.rDa'), package = 'scDNA'))
  
  if(panel=="MSK_RL"){
    print("Loading TxDB for Myeloid Clonal Evolution (Levine, MSK)/MSK_RL")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/MSK_RL_txdb', package = 'scDNA')) 
  } else if(panel=="Myeloid"){
    print("Loading TxDB for Myeloid Panel")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/myeloid_txdb', package = 'scDNA')) 
  } else if(panel=="UCSC"){
    print("Loading TxDB derived from TxDb.Hsapiens.UCSC.hg19.knownGene")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/hg19_ensembl_txdb', package = 'scDNA'))
  } else if(panel=="mm10"){
    print("Loading TxDB derived Mus_musculus.GRCm38.102.gtf")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/mm10_ensembl_txdb', package = 'scDNA')) 
    }
  if(panel=="mm10"){
    load(system.file(paste0('data/mm10_annotation.Rda'), package = 'scDNA'))
    genes_found<-genes(custom_txdb)$gene_id
    complete_gene_annotation<-mm10_annotation%>%
      dplyr::filter(ensemble_geneID%in%genes_found)%>%
      dplyr::select(hgnc_symbol=Symbol,ensembl_canonical_gene=ensemble_geneID,ensemble_txID)%>%
      dplyr::mutate(final_transcript_id=ensemble_txID)
      seq_source<-BSgenome.Mmusculus.UCSC.mm10::Mmusculus

  } else {
  genes_found<-genes(custom_txdb)$gene_id
  complete_gene_annotation<-cBioPortal_annotation%>%
    dplyr::filter(ensembl_canonical_gene%in%genes_found)%>%
    dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
    dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
 seq_source<-BSgenome.Hsapiens.UCSC.hg19::Hsapiens
   }
  print("Load annotation data")
  
  print("Extracting Variant Matrix")
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
  
  #necessary for meaningful GRangess
  SNV_mat$REF<-as(SNV_mat$REF, "DNAStringSet")
  SNV_mat$ALT<-as(SNV_mat$ALT, "DNAStringSet")
  variant_gRange<-GenomicRanges::makeGRangesFromDataFrame(SNV_mat,
                                                          seqnames.field = "CHROM",
                                                          start.field="POS",
                                                          end.field="POS",
                                                          keep.extra.columns=TRUE)
  
  print("Annotating Variants based on location")
  gene_subset<-GenomicFeatures::genes(custom_txdb)
  non_gene_variants<-variant_gRange[-queryHits(findOverlaps(variant_gRange,gene_subset))]
  genic_variant_gRange_subset<-variant_gRange[queryHits(findOverlaps(variant_gRange,gene_subset))]
  
  print(paste("n =",length(non_gene_variants), "variants were not annotated to be in a gene body"))
  print("They can be found in the following regions and will be annotated with genomic location only")
  print(GenomicRanges::reduce(non_gene_variants, min.gapwidth = 50)%>%data.frame())
  
  exon_subset<-GenomicFeatures::exons(custom_txdb)
  exonic_variant_gRange_subset<-genic_variant_gRange_subset[S4Vectors::queryHits(GenomicRanges::findOverlaps(genic_variant_gRange_subset,exon_subset))]
  non_exonic_variant_gRange_subset<-genic_variant_gRange_subset[-S4Vectors::queryHits(GenomicRanges::findOverlaps(genic_variant_gRange_subset,exon_subset))]
  
  print(paste("The following n =",length(genic_variant_gRange_subset), "variants were found within the following regions of a gene body"))
  all_genic_variant_lists<-list(
    "Coding" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::CodingVariants()))),
    "Splice" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::SpliceSiteVariants()))),
    "Intronic" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::IntronVariants()))),
    "FiveUTR" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::FiveUTRVariants()))),
    "ThreeUTR" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::ThreeUTRVariants()))),
    "Promoter" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (genic_variant_gRange_subset), subject = custom_txdb,region = VariantAnnotation::IntergenicVariants()))))
  
  variant_QUERYID_by_region<-lapply(names(all_genic_variant_lists),function(x){
    
    if(length(unique(all_genic_variant_lists[[x]]$QUERYID))>0){
      data.frame("Region"=x,
                 "QUERYID"=(all_genic_variant_lists[[x]]$QUERYID),
                 "LOCSTART"=(all_genic_variant_lists[[x]]$LOCSTART))%>%
        dplyr::distinct(Region,QUERYID,LOCSTART)
      
    } else{
      print(paste("No",x,"variants found"))
      return(NULL)
    }
  })
  
  variant_QUERYID_by_region_df<-do.call(rbind,variant_QUERYID_by_region)
  variant_breakdown<-variant_QUERYID_by_region_df%>%
    dplyr::mutate(Region=factor(Region,levels=names(all_genic_variant_lists)))%>%
    dplyr::group_by(Region)%>%
    dplyr::summarise(Count=dplyr::n())
  variant_prioritized_grouping<-variant_QUERYID_by_region_df%>%
    dplyr::mutate(Region=factor(Region,levels=names(all_genic_variant_lists)))%>%
    dplyr::arrange(Region)%>%
    dplyr::distinct(QUERYID,.keep_all = TRUE)
  variant_breakdown_final<-variant_prioritized_grouping%>%
    dplyr::group_by(Region)%>%
    dplyr::summarise(Count=dplyr::n())%>%
    dplyr::full_join(variant_breakdown,.,by="Region",suffix=c("_Initial","_Reduced"))
  print("Prioritizing annotation for variants that appear in more than one group")
  print(variant_breakdown_final)
  
  
  print("Annotating coding variants")
  variant_annotation_location_list<-setNames(lapply(list("Coding","Splice","Intronic"),function(Region_test){
    region_QUERYID<-variant_prioritized_grouping%>%
      dplyr::filter(Region==Region_test)%>%
      dplyr::pull(QUERYID)
    
    region_gRange<- variant_gRange[variant_gRange$QUERYID%in%region_QUERYID]
    exonic_region_variant_gRanges<-region_gRange[unique(queryHits(findOverlaps(region_gRange,exonic_variant_gRange_subset)))]
    non_exonic_region_variant_gRanges<-region_gRange[unique(queryHits(findOverlaps(region_gRange,non_exonic_variant_gRange_subset)))]
    border_region_ids <- setdiff(region_gRange$id,union(exonic_region_variant_gRanges$id,non_exonic_region_variant_gRanges$id))
    border_region_gRange<-variant_gRange[variant_gRange$id%in%border_region_ids]
    
    exonic_variants <- suppressWarnings(suppressMessages(VariantAnnotation::predictCoding(query = exonic_region_variant_gRanges, 
                                                                                          subject = custom_txdb, 
                                                                                          seqSource = seq_source, 
                                                                                          varAllele = exonic_region_variant_gRanges$ALT)))
    non_exonic_variants <- suppressWarnings(suppressMessages(VariantAnnotation::predictCoding(query = non_exonic_region_variant_gRanges, 
                                                                                              subject = custom_txdb, 
                                                                                              seqSource = seq_source, 
                                                                                              varAllele = non_exonic_region_variant_gRanges$ALT)))
    border_variants <- suppressWarnings(suppressMessages(VariantAnnotation::predictCoding(query = border_region_gRange, 
                                                                                          subject = custom_txdb, 
                                                                                          seqSource = seq_source, 
                                                                                          varAllele = border_region_gRange$ALT)))
    
    list("Exonic"=exonic_variants,
         "Non_exonic"=non_exonic_variants,
         "Border"=border_variants)
  }),c("Coding","Splice","Intronic"))
  
  print("Variants successfully annotated for impact on coding change")
  print("Only a subset of coding variants will have a final impact on protein sequence")
  print(do.call(rbind,lapply(variant_annotation_location_list,function(type){
    lapply(type,function(region){
      length(region)
    })
  })))
  
  final_protein_annotation<-do.call(rbind,lapply(unlist(variant_annotation_location_list)%>%names,function(variant_list){
    unlist(variant_annotation_location_list)[[variant_list]]%>%
      data.frame%>%
    dplyr::mutate(Class=variant_list)}))%>%
    dplyr::full_join(complete_gene_annotation,by=c("GENEID"="ensembl_canonical_gene"))%>%
    dplyr::mutate(AA_change=paste0(hgnc_symbol,".",REFAA,PROTEINLOC,VARAA))#%>%
  
  non_annotated_genic_GRanges<-genic_variant_gRange_subset[!genic_variant_gRange_subset$id%in%final_protein_annotation$id]
  non_annotated_genic_GRanges$GENEID<-gene_subset[findOverlaps(non_annotated_genic_GRanges,gene_subset,select="first")]$gene_id
  non_annotated_nongenic_GRanges<-genic_variant_gRange_subset[!genic_variant_gRange_subset$id%in%final_protein_annotation$id]
  
  
  final_annotation<- final_protein_annotation%>%
    dplyr::bind_rows(non_annotated_genic_GRanges%>%
                       data.frame%>%
                       dplyr::mutate(Class="non_coding")%>%
                       dplyr::inner_join(complete_gene_annotation,by=c("GENEID"="ensembl_canonical_gene")))%>%
    dplyr::filter(!is.na(id))%>%
    dplyr::bind_rows(non_gene_variants%>%
                       data.frame%>%
                       dplyr::mutate(Class="Intergenic"))%>%
    dplyr::full_join(variant_prioritized_grouping%>%
                       dplyr::select(QUERYID,LOCSTART),
                     by="QUERYID")%>%
    dplyr::select(id,seqnames,start,end,width,strand,REF,ALT,QUAL,amplicon,QUERYID,
                  final_transcript_id,SYMBOL=hgnc_symbol,CDSLOC=LOCSTART, REFAA,PROTEINLOC,VARAA,AA_change,CONSEQUENCE,Class)%>%
    dplyr::mutate(CDS_change=paste0(SYMBOL,":c.",CDSLOC,REF,">",ALT))%>%
    dplyr::mutate(final_annot=dplyr::case_when(
      !is.na(AA_change)~AA_change,
      is.na(AA_change)&!is.na(CDS_change)~CDS_change,
      is.na(AA_change)&is.na(CDS_change)~id
    ))%>%
    dplyr::mutate(final_annot=dplyr::case_when(
      !is.na(AA_change)~AA_change,
      !is.na(SYMBOL)&!is.na(REFAA)~CDS_change,
      #     is.na(AA_change)&!is.na(CDS_change)~CDS_change,
      Class=="non_coding"~CDS_change,
      is.na(SYMBOL)~id,
      TRUE~id
    ))%>%
    dplyr::mutate(Class=dplyr::case_when(
      Class=="Coding.Exonic"~"Exon",
      Class=="Intronic.Exonic"~"Exon_Boundary",
      Class=="non_coding"~"Intronic",
      TRUE~Class
    ))
  
  print("Final annotation complete")
  return(final_annotation)
}
