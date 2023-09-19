#' Title
#'
#' @param panel_bed_path 
#'
#' @return
#' @export
#'
#' @examples
gemerate_txdb<-function(panel_bed_path,
                        panel_name=NULL,
                        save_path=NULL){
    
    print("Extracting all genes from BED file")
    genes<-genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    bed<-read.delim(panel_bed_path,sep="\t",header=FALSE)%>%
      dplyr::select(chr=V1,start=V2,end=V3,amplicon=V4)%>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
    
    genes_found<-genes[subjectHits(findOverlaps(bed,genes))]$gene_id
    load(system.file('data/cBioPortal_annotation.rDa', package = 'scDNA')) # loads in as variable annotation_file?
    
    complete_gene_annotation<-cBioPortal_annotation%>%
                dplyr::filter(entrez_gene_id%in%genes_found)%>%
                dplyr::select(hgnc_symbol,ensembl_canonical_gene,ensembl_canonical_transcript,mskcc_canonical_transcript,ccds_id)%>%
                dplyr::mutate(final_transcript_id=ifelse(mskcc_canonical_transcript=="",ensembl_canonical_transcript,mskcc_canonical_transcript))
              
    print(paste0("Identified n = ",length(unique(complete_gene_annotation$hgnc_symbol))," genes:"))
    print(unique(complete_gene_annotation$hgnc_symbol))
          
    print("Generating custom TxDB")
    custom_txdb<-GenomicFeatures::makeTxDbFromUCSC(genome="hg19", tablename="ensGene",
                                                 transcript_ids=complete_gene_annotation$final_transcript_id,
                                                 circ_seqs=NULL,
                                                 url="http://genome.ucsc.edu/cgi-bin/",
                                                 goldenPath.url=getOption("UCSC.goldenPath.url"),
                                                 taxonomyId=NA,
                                                 miRBaseBuild=NA)
    saveDb(custom_txdb,file=save_path)
    return(custom_txdb)
}