
#' Normalize Protein Data
#'
#' @param file
#' @param metadata
#' @param protein_mat
#' @param method
#' @param detect_IgG
#' @param NGT
#'
#' @return
#' @export
#'
#' @examples
normalize_protein_data<-function(file,
                                 metadata,
                                 protein_mat,
                                 method="dsb",
                                 detect_IgG=TRUE,
                                 background_droplets,
                                 NGT){

  if(method=="dsb"){

    all_protein_droplets <- rhdf5::h5read(file=file,
                                          name="/all_barcodes/protein_read_counts/layers/read_counts")
    all_dna_droplets <- rhdf5::h5read(file=file,
                                      name="/all_barcodes/dna_read_counts/layers/read_counts")
    colnames(all_dna_droplets) <-rhdf5::h5read(file=file,
                                               name="/all_barcodes/dna_read_counts/ra/barcode")
    colnames(all_protein_droplets) <- rhdf5::h5read(file=file,
                                                    name="/all_barcodes/protein_read_counts/ra/barcode")
    colnames(all_protein_droplets) <- gsub("-1","",colnames(all_protein_droplets))


    cell_protein_matrix_input <- t(data.frame(protein_mat%>%
                                                dplyr::select(!Cell)))

    empty_drops_matrix_input <- data.frame(all_protein_droplets)%>%
                                    dplyr::select(all_of(background_droplets))

    isotype <- grep("IgG",colnames(protein_mat),value=TRUE)


    adt_norm <- DSBNormalizeProtein(
      # remove ambient protien noise reflected in counts from empty droplets
      cell_protein_matrix = cell_protein_matrix_input, # cell-containing droplet raw protein count matrix
      empty_drop_matrix = empty_drops_matrix_input, # empty/background droplet raw protein counts

      # recommended step II: model and remove the technical component of each cell's protein library
      denoise.counts = TRUE, # (default = TRUE); run step II
      use.isotype.control = TRUE, # (default = TRUE): use isotype controls to define technical components.
      isotype.control.name.vec = isotype#,# vector of isotype control names
    )
  } else{
    print("method not available")
    break
  }
}
