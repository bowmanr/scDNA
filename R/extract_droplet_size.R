

#' Extract protein library size, dna library size, and amplicon size for all droplets.
#'
#' @param file path to h5 file
#' @param sce singleCellExperiment Object
#' @return metadata matrix
#' @export
#'
#' @examples
extract_droplet_size<- function(sce){
  file<-sce@metadata$file
  all_protein_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/layers/read_counts")
  all_dna_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/layers/read_counts")
  colnames(all_dna_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/ra/barcode")
  colnames(all_protein_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/ra/barcode")

  # create a metadata dataframe of simple qc stats for each droplet
  dna_size <- data.frame("Cell"=colnames(all_dna_droplets),
                        "dna_size"=log10(Matrix::colSums(all_dna_droplets)),
                        "amplicons"=Matrix::colSums(all_dna_droplets > 0))
  protein_size <- data.frame("Cell"=colnames(all_protein_droplets),
                             "proteins"=Matrix::colSums(all_protein_droplets > 0),
                            "protein_size"=log10(Matrix::colSums(all_protein_droplets)))
  md <- dplyr::inner_join(dna_size,protein_size)%>%
         dplyr::mutate(Cell=gsub("-1","",Cell))
  md<- md%>%
          mutate(Droplet_type=case_when(
                    Cell%in%colnames(sce)~"Cell",
                    TRUE~"Empty"
                  ))%>%
          full_join(sce@metadata$NGT,by="Cell")
  return(md)
}
