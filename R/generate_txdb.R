#' generating a transcription data base (TxDb)
#' @description
#' This helper function allows for custom TxDbs to be generated. Note: our annotation method relies on cBioPortal's transcripts as a recommended transcript of interest.
#'
#' @param gtf_gz_file_path a gtf.gz file for creating a custom txDb and suggested cbioPortalTranscript dataframe.
#' @param folder_saved_directory where to save the custom TxDB file and suggested cbioportal_files.
#'
#' @return It returns a list of custom TxDb file and a recommended cbioPortal.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' custom_txdb <- generate_txdb(panel_bed_path="~/my_bed_file.bed",folder_saved_directory= "~/my_custom_txdb")
#' }
generate_txdb<-function(gtf_gz_file_path,
                        folder_saved_directory=NULL){


  gtf <- data.table::fread(
    gtf_gz_file_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    col.names = c(
      "seqname", "source", "feature", "start", "end",
      "score", "strand", "frame", "attribute"
    )
  )
  strip_gtf_attribute_version <- function(attribute, key) {
    pattern <- paste0('(', key, ' ")([^"]+?)\\.[0-9]+(")')
    replacement <- "\\1\\2\\3"
    stringr::str_replace_all(attribute, pattern, replacement)
  }
  gtf <- gtf%>%
    dplyr::mutate(
      attribute = strip_gtf_attribute_version(attribute, "gene_id"),
      attribute = strip_gtf_attribute_version(attribute, "transcript_id"),
      attribute = strip_gtf_attribute_version(attribute, "exon_id"),
    )

  custom_gtf_file <- paste0(folder_saved_directory,"custom_gtf.gtf")

  data.table::fwrite(
    gtf %>%
      dplyr::select(
        seqname, source, feature, start, end,
        score, strand, frame, attribute
      ),
    file = custom_gtf_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )

  custom_txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = custom_gtf_file,
    format = "gtf",
    dataSource = "From Custom GTF",
    organism = "Homo sapiens",
    taxonomyId = 9606,
    circ_seqs = "chrM"
  )

  GenomeInfoDb::seqlevelsStyle(custom_txdb) <- "UCSC"
  AnnotationDbi::saveDb(custom_txdb,paste0(folder_saved_directory,"hg38_custom_txdb"))

  gtf <- gtf %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::mutate(
      gene_id = stringr::str_match(attribute, 'gene_id "([^"]+)"')[,2],
      transcript_id =  stringr::str_match(attribute, 'transcript_id "([^"]+)"')[,2],
      exon_id =  stringr::str_match(attribute, 'exon_id "([^"]+)"')[,2],
      gene_name =  stringr::str_match(attribute, 'gene_name "([^"]+)"')[,2],
      gene_type =  stringr::str_match(attribute, 'gene_type "([^"]+)"')[,2],
      transcript_type =  stringr::str_match(attribute, 'transcript_type "([^"]+)"')[,2],
      transcript_name =  stringr::str_match(attribute, 'transcript_name "([^"]+)"')[,2],
      tag_all = attribute,
      gene_id_base = sub("\\..*$", "", gene_id),
      transcript_id_base = sub("\\..*$", "", transcript_id),
      exon_width = end - start + 1,
      is_mane_select =  stringr::str_detect(tag_all, 'tag "MANE_Select"'),
      is_mane_plus_clinical =  stringr::str_detect(tag_all, 'tag "MANE_Plus_Clinical"'),
      is_ensembl_canonical =  stringr::str_detect(tag_all, 'tag "Ensembl_canonical"'),
      is_gencode_basic =  stringr::str_detect(tag_all, 'tag "basic"'),
      is_appris_principal =  stringr::str_detect(tag_all, 'tag "appris_principal'))

  load(system.file(paste0('data/cBioPortal_annotation.rDa'), package = 'scDNA'))

  cBioPortal_annotation<-cBioPortal_annotation%>%
    dplyr::left_join(gtf%>%
                dplyr::select(seqname,source,feature,start,end,exon_width,gene_name,transcript_name,transcript_type,gene_id_base,transcript_id_base,is_ensembl_canonical,is_mane_select,score, strand, frame, attribute ), by=c("hgnc_symbol"="gene_name"))%>%
    dplyr::select(seqname,source,feature,start,end,hgnc_symbol,transcript_name,transcript_type,gene_id_base,ensembl_canonical_gene,transcript_id_base,mskcc_canonical_transcript,is_ensembl_canonical,score, strand, frame, attribute,exon_width)%>%
    dplyr::group_by(transcript_id_base)%>%
    dplyr::mutate(total_exon_width = sum(exon_width))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(hgnc_symbol)%>%
    dplyr::mutate(longest_tx_flag = total_exon_width == max(total_exon_width),
           match_flag = mskcc_canonical_transcript == transcript_id_base,
           ensembl_flag = is_ensembl_canonical)%>%
    dplyr::distinct(hgnc_symbol,transcript_name,gene_id_base,ensembl_canonical_gene,transcript_id_base,mskcc_canonical_transcript,longest_tx_flag,match_flag,ensembl_flag)%>%
    dplyr::arrange( dplyr::desc(match_flag), dplyr::desc(longest_tx_flag), dplyr::desc(ensembl_flag))%>%
    dplyr::mutate(rank_to_keep =  dplyr::row_number())%>%
    dplyr::mutate(suggested = ifelse(rank_to_keep==1,TRUE,FALSE))%>%
    dplyr::ungroup()%>%
    dplyr::select(hgnc_symbol,transcript_name,gene_id_base,ensembl_canonical_gene,transcript_id_base,mskcc_canonical_transcript,longest_tx_flag,match_flag,suggested)

  save(cBioPortal_annotation,file = paste0(folder_saved_directory,"cBioportal_custom_annotation.Rda"))

  return(list("custom_txdb" =custom_txdb,
              "recommended_transcripts"=cBioPortal_annotation))

  }
