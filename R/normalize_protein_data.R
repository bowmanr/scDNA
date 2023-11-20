
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
normalize_protein_data<-function(sce,
                                 metadata,
                                 method="dsb",
                                 detect_IgG=TRUE,
                                 background_droplets){
  file<-sce@metadata$file
  protein_sce <- SingleCellExperiment::altExp(sce,"Protein")
  protein_mat <- protein_sce@assays@data$Protein
  if("DSB"%in%method|"dsb"%in%method){
    print("DSB normalization")
    cells_of_interest<-colnames(protein_mat)
    
    all_protein_droplets <- rhdf5::h5read(file=file,
                                          name="/all_barcodes/protein_read_counts/layers/read_counts")
    colnames(all_protein_droplets) <- rhdf5::h5read(file=file,
                                          name="/all_barcodes/protein_read_counts/ra/barcode")
    colnames(all_protein_droplets) <- gsub("-1","",colnames(all_protein_droplets))
    

    colnames(all_protein_droplets) <-gsub("-1","",colnames(all_protein_droplets))
    
    empty_drops_matrix_input <- data.frame(all_protein_droplets) %>% 
                                    dplyr::select(all_of(background_droplets))
    
    rownames(empty_drops_matrix_input)<- rownames(protein_mat)
    
    if(detect_IgG){
        isotype <- grep("IgG",colnames(protein_mat),value=TRUE)
    }else{
        isotype <- FALSE
    }

    adt_norm <- dsb::DSBNormalizeProtein(
      # remove ambient protien noise reflected in counts from empty droplets
      cell_protein_matrix = protein_mat, # cell-containing droplet raw protein count matrix
      empty_drop_matrix = empty_drops_matrix_input, # empty/background droplet raw protein counts
      # recommended step II: model and remove the technical component of each cell's protein library
      denoise.counts = TRUE, # (default = TRUE); run step II
      use.isotype.control = detect_IgG, # (default = TRUE): use isotype controls to define technical components.
      isotype.control.name.vec = isotype#,# vector of isotype control names
    )
    SingleCellExperiment::assay(protein_sce, "DSB_norm")<-adt_norm
    
  } 
  if("CLR"%in%method|"clr"%in%method){
    print("CLR normalization")
    s <- Seurat::CreateSeuratObject(counts=protein_mat%>%Seurat::as.sparse(), 
                            assay="Protein")
    s <- Seurat::NormalizeData(s,normalization.method = "CLR")
    SingleCellExperiment::assay(protein_sce, "CLR_norm")<-s@assays$Protein@data
    
  }
  
  if(!all(method%in%c("dsb","CLR")) ){
    print("method not available")
    break
  }
  protein_sce@metadata<-metadata
  protein_sce@colData<-S4Vectors::DataFrame(metadata%>%
                                            dplyr::filter(Cell%in%colnames(protein_mat)))
  SingleCellExperiment::altExp(sce, "Protein") <- protein_sce
  
  return(sce)
}
