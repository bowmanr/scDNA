#' Extract Droplet Size
#' @description
#' Extract protein library size, dna library size, and amplicon size for all droplets.
#'
#' @param sce singleCellExperiment object
#'
#' @return metadata dataframe which contains droplet_name, dna_size, amplicon, protein_id, and protein_size as well as the NGT matrix
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' droplet_metadata<- extract_droplet_size(sce)
#' }
extract_droplet_size<- function(sce){
  file<-sce@metadata$file
  all_protein_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/layers/read_counts")
  all_dna_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/layers/read_counts")
  colnames(all_dna_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/ra/barcode")
  colnames(all_protein_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/ra/barcode")

  dna_size <- data.frame("Cell"=colnames(all_dna_droplets),
                        "dna_size"=log10(Matrix::colSums(all_dna_droplets)),
                        "amplicons"=Matrix::colSums(all_dna_droplets > 0))
  protein_size <- data.frame("Cell"=colnames(all_protein_droplets),
                             "proteins"=Matrix::colSums(all_protein_droplets > 0),
                            "protein_size"=log10(Matrix::colSums(all_protein_droplets)))
  md <- dplyr::inner_join(dna_size,protein_size)%>%
         dplyr::mutate(Cell=gsub("-1","",Cell))%>%
    dplyr::mutate(Droplet_type=dplyr::case_when(
           Cell%in%colnames(sce)~"Cell",
           TRUE~"Empty"
         ))

  if("NGT"%in%names(sce@metadata)){
    message("NGT found in sce. Adding to metadata")
    md<- md%>%dplyr::full_join(sce@metadata$NGT,by="Cell")
  }

  return(md)
}
