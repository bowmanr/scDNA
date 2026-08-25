#' Annotate variants of interest
#' @description
#' This function takes in h5 files to extract DNA information through given TXDB files (primarily hg19)
#' The gene names are mapped for their specific nucleotide positions. The starting and ending positions
#' for the chromosome are listed. The consequence of the mutation such as synonymous, nonsynonymous,
#' as well as if this is a coding region, on the exon boundry or intronic. Short amino acid changes are also labeled.
#' The variance matrix also extracts the amplicon that the variant is found on.
#'
#' @param file path to h5 file to pull out relevant DNA information
#' @param panel Which panel to use often passed from variant_ID function
#' @param select_variants variants of interest, default is all variants
#'
#' @importFrom methods as
#' @importFrom rlang .data
#' @return variant annotation matrix
#' @export
#' @examples
#' \dontrun{
#' variant_table <- annotate_variants(file= "./Data/sample_file.h5",panel="hg19")
#' }
annotate_variants<- function(file,
                             panel=NULL,
                             select_variants=NULL){

  print("Loading annotation data")
  if(panel=="MSK_RL"){
    print("Loading TxDB for Myeloid Clonal Evolution (Levine, MSK)/MSK_RL")
    load(system.file(paste0('data/cBioPortal_annotation.rDa'), package = 'scDNA'))
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/MSK_RL_txdb', package = 'scDNA'))
    complete_gene_annotation<-cBioPortal_annotation%>%
      dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
      dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))%>%
      dplyr::mutate(transcript_id_base=final_transcript_id)%>%
      dplyr::mutate(gene_id_base=ensembl_canonical_gene)

    seq_source<-BSgenome.Hsapiens.UCSC.hg19::Hsapiens

  } else if(panel=="Myeloid"){
    print("Loading TxDB for Myeloid Panel")
    load(system.file(paste0('data/cBioPortal_annotation.rDa'), package = 'scDNA'))
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/myeloid_txdb', package = 'scDNA'))
    complete_gene_annotation<-cBioPortal_annotation%>%
      dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
      dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))%>%
    dplyr::mutate(transcript_id_base=final_transcript_id)%>%
    dplyr::mutate(gene_id_base=ensembl_canonical_gene)
    seq_source<-BSgenome.Hsapiens.UCSC.hg19::Hsapiens

  } else if(panel=="UCSC"){
    print("Loading TxDB derived from gencode.v19.annotation.gtf.gz")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/hg19_custom_txdb', package = 'scDNA'))
    load(system.file(paste0('data/cBioPortal_hg19_annotation.rDa'), package = 'scDNA'))
    complete_gene_annotation<-cBioPortal_annotation_hg19%>%
      dplyr::select(hgnc_symbol,gene_id_base,mskcc_canonical_transcript,transcript_id_base,suggested)%>%
      dplyr::mutate(final_transcript_id=transcript_id_base )
    seq_source<-BSgenome.Hsapiens.UCSC.hg19::Hsapiens

  } else if(panel=="mm10"){
    print("Loading TxDB derived Mus_musculus.GRCm38.102.gtf")
    custom_txdb<- AnnotationDbi::loadDb(system.file('data/mm10_ensembl_txdb', package = 'scDNA'))
    load(system.file(paste0('data/mm10_annotation.Rda'), package = 'scDNA'))
    genes_found<-GenomicFeatures::genes(custom_txdb)$gene_id
    complete_gene_annotation<-mm10_annotation%>%
      dplyr::select(hgnc_symbol=Symbol,ensembl_canonical_gene=ensemble_geneID,ensemble_txID)%>%
      dplyr::mutate(final_transcript_id=ensemble_txID)%>%
      dplyr::mutate(transcript_id_base=final_transcript_id)


    seq_source<-BSgenome.Mmusculus.UCSC.mm10::Mmusculus

  } else if(panel=="hg38"){
  print("Loading TxDB derived from gencode.v50.annotation.gtf.gz")
  custom_txdb<- AnnotationDbi::loadDb(system.file('data/hg38_custom_txdb', package = 'scDNA'))

  GenomeInfoDb::genome(custom_txdb)<-"hg38"
  load(system.file(paste0('data/cBioportal_hg38_annotation.Rda'), package = 'scDNA'))
  complete_gene_annotation<-cBioPortal_annotation_hg38%>%
    dplyr::select(hgnc_symbol,gene_id_base,mskcc_canonical_transcript,transcript_id_base,suggested)%>%
    dplyr::mutate(final_transcript_id=transcript_id_base )
  seq_source <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  }


  print("Extracting Variant Matrix")
  if(grepl("loom",file)){
    SNV_mat_prefilter<-data.frame(do.call(cbind, rhdf5::h5read(file=file,name="/row_attrs/")))
  } else if(grepl("h5",file)){
    SNV_mat_prefilter<-data.frame(do.call(cbind, rhdf5::h5read(file=file,name="/assays/dna_variants/ca/")))
  }

  SNV_mat_prefilter<- SNV_mat_prefilter%>%
    dplyr::select(id,ALT,REF,CHROM,POS,QUAL,amplicon)%>%
    dplyr::mutate(ALT=gsub("\\*","N",.data$ALT))%>%
    dplyr::mutate(REF=gsub("\\*","N",.data$REF))%>%
    dplyr::mutate(CHROM=paste0("chr",.data$CHROM))%>%
    dplyr::mutate(QUAL=as.numeric(QUAL))%>%
    dplyr::mutate(QUERYID=1:nrow(.))%>%
    dplyr::mutate(POS=as.numeric(POS))%>%
    dplyr::mutate(width = pmax(nchar(as.character(REF)), 1),
                  END = POS + width - 1  )%>%
    dplyr::mutate(POS=as.character(POS),
                  END=as.character(END))

  SNV_mat <- SNV_mat_prefilter
  
  if (!is.null(select_variants)) {
    SNV_mat <- SNV_mat %>% dplyr::filter(id %in% select_variants)
  }

  SNV_mat$REF<-as(SNV_mat$REF, "DNAStringSet")
  SNV_mat$ALT<-as(SNV_mat$ALT, "DNAStringSet")

  variant_gRange<-GenomicRanges::makeGRangesFromDataFrame(SNV_mat,
                                                          seqnames.field = "CHROM",
                                                          start.field="POS",
                                                          end.field="END",
                                                          keep.extra.columns=TRUE)

  names(variant_gRange) <- variant_gRange$id
  GenomeInfoDb::seqlevelsStyle(variant_gRange) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(custom_txdb) <- "UCSC"

  common_seqlevels <- intersect(GenomeInfoDb::seqlevels(variant_gRange),GenomeInfoDb::seqlevels(custom_txdb))
  variant_gRange <- GenomeInfoDb::keepSeqlevels(variant_gRange,common_seqlevels,pruning.mode = "coarse")
  GenomeInfoDb::seqinfo(variant_gRange) <- GenomeInfoDb::seqinfo(custom_txdb)[common_seqlevels]
  GenomeInfoDb::genome(variant_gRange) <- GenomeInfoDb::genome(custom_txdb)


  print("Annotating Variants based on location")
  gene_subset<-GenomicFeatures::genes(custom_txdb)
  non_gene_variants<-variant_gRange[-S4Vectors::queryHits(GenomicRanges::findOverlaps(variant_gRange,gene_subset))]
  genic_variant_gRange_subset<-variant_gRange[S4Vectors::queryHits(GenomicRanges::findOverlaps(variant_gRange,gene_subset))]

  non_gene_variants<-subsetByOverlaps(variant_gRange,gene_subset,invert = T)
  genic_variant_gRange_subset<-subsetByOverlaps(variant_gRange,gene_subset,invert = F)

  print(paste("n =",length(non_gene_variants), "variants were not annotated to be in a gene body"))
  print("They can be found in the following regions and will be annotated with genomic location only")
  print(GenomicRanges::reduce(non_gene_variants, min.gapwidth = 50)%>%data.frame())

  exon_subset<-GenomicFeatures::exons(custom_txdb)
  exonic_variant_gRange_subset<-variant_gRange[S4Vectors::queryHits(GenomicRanges::findOverlaps(variant_gRange,exon_subset))]
  non_exonic_variant_gRange_subset<-variant_gRange[-S4Vectors::queryHits(GenomicRanges::findOverlaps(variant_gRange,exon_subset))]

  print(paste("The following n =",length(genic_variant_gRange_subset), "variants were found within the following regions of a gene body"))
  all_genic_variant_lists<-list(
    "Coding" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::CodingVariants()))),
    "Splice" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::SpliceSiteVariants()))),
    "Intronic" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::IntronVariants()))),
    "FiveUTR" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::FiveUTRVariants()))),
    "ThreeUTR" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::ThreeUTRVariants()))),
    "Promoter" = suppressWarnings(suppressMessages(VariantAnnotation::locateVariants(query = (variant_gRange), subject = custom_txdb,region = VariantAnnotation::IntergenicVariants()))))


  transcript_subset <- GenomicFeatures::transcripts(custom_txdb)
  tx_lookup <- AnnotationDbi::select(custom_txdb,
                                     keys = as.character(transcript_subset$tx_id),
                                 keytype = "TXID",
                                     columns = c("TXID", "TXNAME", "GENEID")) %>%
    dplyr::distinct() %>%
    dplyr::rename(TXID = TXID,
                  tx_name = TXNAME,
                  GENEID_txdb = GENEID) %>%
    dplyr::mutate(TXID = as.integer(TXID))

  variant_QUERYID_by_region<-lapply(names(all_genic_variant_lists),function(x){

    if(length(unique(all_genic_variant_lists[[x]]$QUERYID))>0){
      data.frame("Region"=x,
                 "variant_id"=names(all_genic_variant_lists[[x]]),
                 "QUERYID"=(all_genic_variant_lists[[x]]$QUERYID),
                 "TXID" = as.numeric(all_genic_variant_lists[[x]]$TXID),
                 "GENEID" = (all_genic_variant_lists[[x]]$GENEID),
                 "LOCSTART"=(all_genic_variant_lists[[x]]$LOCSTART),
                 "LOCEND"=(all_genic_variant_lists[[x]]$LOCEND))%>%
        dplyr::left_join(tx_lookup, by = "TXID") %>%
        dplyr::distinct(Region,variant_id,QUERYID,TXID,tx_name,GENEID,GENEID_txdb,LOCSTART,LOCEND)

    } else{
      print(paste("No",x,"variants found"))
      return(NULL)
    }
  })
  variant_QUERYID_by_region_df<-do.call(rbind,variant_QUERYID_by_region)

  variant_tx_summary <- variant_QUERYID_by_region_df %>%
    dplyr::group_by(QUERYID) %>%
    dplyr::summarise(tx_name = paste(unique(stats::na.omit(tx_name)), collapse = ";"),
                     TXID_list = paste(unique(stats::na.omit(TXID)), collapse = ";"),
                     GENEID_txdb_list = paste(unique(stats::na.omit(GENEID_txdb)), collapse = ";"),
                     .groups = "drop") %>%
    dplyr::mutate(tx_name = dplyr::na_if(tx_name, ""),
                  TXID_list = dplyr::na_if(TXID_list, ""),
                  GENEID_txdb_list = dplyr::na_if(GENEID_txdb_list, ""))

  variant_prioritized_grouping <- variant_QUERYID_by_region_df %>%
    dplyr::mutate(Region = factor(Region, levels = names(all_genic_variant_lists))) %>%
    dplyr::arrange(Region) %>%
    dplyr::distinct(QUERYID, .keep_all = TRUE) %>%
    dplyr::select(Region,variant_id,QUERYID,LOCSTART,LOCEND)%>%
    dplyr::left_join(variant_tx_summary, by = "QUERYID")



  variant_breakdown<-variant_QUERYID_by_region_df%>%
    dplyr::mutate(Region=factor(Region,levels=names(all_genic_variant_lists)))%>%
    dplyr::group_by(Region)%>%
    dplyr::summarise(Count=dplyr::n())
  variant_breakdown_final<-variant_prioritized_grouping%>%
    dplyr::group_by(Region)%>%
    dplyr::summarise(Count=dplyr::n())%>%
    dplyr::full_join(variant_breakdown,.,by="Region",suffix=c("_Initial","_Reduced"))
  print("Prioritizing annotation for variants that appear in more than one group")
  print(variant_breakdown_final)


  print("Annotating coding variants")
  variant_annotation_location_list<-setNames(lapply(list("Coding","Splice","Intronic","Promoter"),function(Region_test){
  region_QUERYID<-variant_QUERYID_by_region_df%>%
      dplyr::filter(Region==Region_test)%>%
      dplyr::pull(variant_id)

    region_indices <- match(region_QUERYID, names(variant_gRange))
    region_gRange <- variant_gRange[region_indices[!is.na(region_indices)]]

    exonic_region_variant_gRanges<-region_gRange[unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(region_gRange,exonic_variant_gRange_subset)))]
    non_exonic_region_variant_gRanges<-region_gRange[unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(region_gRange,non_exonic_variant_gRange_subset)))]


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
  }),c("Coding","Splice","Intronic","Promoter"))

  print("Variants successfully annotated for impact on coding change")
  print("Only a subset of coding variants will have a final impact on protein sequence")
  print(do.call(rbind,lapply(variant_annotation_location_list,function(type){
    lapply(type,function(region){
      length(region)
    })
  })))

    final_protein_annotation <- do.call(rbind,lapply(names(unlist(variant_annotation_location_list)),function(variant_list) {
          x <- unlist(variant_annotation_location_list)[[variant_list]]
          if (length(x) == 0) {
            return(NULL)
          }
          x %>%
            data.frame() %>%
            dplyr::mutate(Class = variant_list,
                          TXID=as.integer(TXID))}))%>%
      dplyr::left_join(tx_lookup,by = "TXID") %>%
      dplyr::left_join(complete_gene_annotation,by = c("tx_name" = "transcript_id_base")) %>%
      dplyr::mutate(tx_name = dplyr::coalesce(tx_name, final_transcript_id),
                    SYMBOL = hgnc_symbol,
                    AA_change = dplyr::case_when(
                      !is.na(SYMBOL) & !is.na(REFAA) & !is.na(PROTEINLOC) & !is.na(VARAA) ~ paste0(SYMBOL, ".", REFAA, PROTEINLOC, VARAA),
                      TRUE ~ NA_character_))


    non_annotated_genic_GRanges <- genic_variant_gRange_subset[!names(genic_variant_gRange_subset) %in% unique(final_protein_annotation$id)]
    gene_hits <- GenomicRanges::findOverlaps(non_annotated_genic_GRanges,gene_subset,select = "first")
    non_annotated_genic_GRanges$GENEID <- gene_subset[gene_hits]$gene_id

    tx_hits <- GenomicRanges::findOverlaps(non_annotated_genic_GRanges,transcript_subset)

    noncoding_tx_summary <- data.frame(id = non_annotated_genic_GRanges$id[S4Vectors::queryHits(tx_hits)],
                                       TXID = transcript_subset$tx_id[S4Vectors::subjectHits(tx_hits)])%>%
      dplyr::left_join(tx_lookup, by = "TXID") %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        tx_name = paste(unique(stats::na.omit(tx_name)), collapse = ";"),
        TXID_list = paste(unique(stats::na.omit(TXID)), collapse = ";"),
        GENEID_txdb_list = paste(unique(stats::na.omit(GENEID_txdb)), collapse = ";"),
        .groups = "drop") %>%
      dplyr::mutate(tx_name = dplyr::na_if(tx_name, ""),
                    TXID_list = dplyr::na_if(TXID_list, ""),
                    GENEID_txdb_list = dplyr::na_if(GENEID_txdb_list, ""))

    noncoding_annotation <- non_annotated_genic_GRanges%>%
      data.frame()%>%
      dplyr::mutate(Class = "non_coding")%>%
      dplyr::left_join(noncoding_tx_summary,by = "id")%>%
      dplyr::left_join(complete_gene_annotation,by = c("tx_name" = "transcript_id_base")) %>%
      dplyr::mutate(SYMBOL = hgnc_symbol,
                    final_transcript_id = dplyr::coalesce(tx_name, final_transcript_id))



    final_annotation <- final_protein_annotation %>%
      dplyr::bind_rows(noncoding_annotation) %>%
      dplyr::filter(!is.na(id)) %>%
      dplyr::bind_rows(non_gene_variants %>%
                         data.frame() %>%
          dplyr::mutate(Class = "Intergenic",
                        tx_name = NA_character_,
                        SYMBOL = NA_character_,
                        final_transcript_id = NA_character_,
                        AA_change = NA_character_,
                        CONSEQUENCE = NA_character_,
                        QUERYID = as.integer(NA),
                        TXID_list = NA_character_,
                        GENEID_txdb_list = NA_character_)) %>%

      dplyr::left_join(variant_prioritized_grouping %>%
                         dplyr::select(variant_id,QUERYID_prioritized = QUERYID,Region,LOCSTART,LOCEND,tx_name_prioritized = tx_name,TXID_list_prioritized = TXID_list,GENEID_txdb_list_prioritized = GENEID_txdb_list),
                       by = c("id" = "variant_id")) %>%
      dplyr::mutate(QUERYID = dplyr::coalesce(as.integer(QUERYID),as.integer(QUERYID_prioritized)),
                    tx_name = dplyr::coalesce(tx_name,tx_name_prioritized),
                    TXID_list = dplyr::coalesce(as.character(TXID_list),as.character(TXID_list_prioritized)),
                    GENEID_txdb_list = dplyr::coalesce(as.character(GENEID_txdb_list),as.character(GENEID_txdb_list_prioritized)),
                    final_transcript_id = dplyr::coalesce(final_transcript_id,tx_name),
                    CDSLOC = LOCSTART,
                    CDS_change = dplyr::case_when(
                      !is.na(SYMBOL) & !is.na(CDSLOC) ~ paste0(SYMBOL, ":c.", CDSLOC, REF, ">", ALT),
                      TRUE ~ NA_character_),
                    final_annot = dplyr::case_when(
                      !is.na(AA_change) ~ AA_change,
                      !is.na(SYMBOL) & !is.na(CDS_change) ~ CDS_change,
                      !is.na(tx_name) ~ paste0(id, ":", tx_name),
                      TRUE ~ id),
                    Class = dplyr::case_when(
                      Class == "Coding.Exonic" ~ "Exon",
                      Class == "Intronic.Exonic" ~ "Exon_Boundary",
                      Class == "non_coding" ~ "Intronic",
                      TRUE ~ Class)) %>%
      dplyr::select(-QUERYID_prioritized,-tx_name_prioritized,-TXID_list_prioritized,-GENEID_txdb_list_prioritized,-QUERYID.1,-QUERYID, -TXID, -CDSID, -GENEID_txdb, -gene_id_base, -tx_name,-TXID_list, -GENEID_txdb_list)%>%
      dplyr::distinct()

  print("Final annotation complete")
  return(final_annotation)

}
