
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
                                 num_components_to_keep=NULL,
                                 background_droplets){
  protein_sce <- SingleCellExperiment::altExp(sce,"Protein")
  protein_mat <- protein_sce@assays@data$Protein
  if("DSB"%in%method|"dsb"%in%method){
    print("DSB normalization")
    file<-sce@metadata$file
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
    SummarizedExperiment::assay(protein_sce, "DSB_norm")<-adt_norm
    
  } 
  if ("CLR" %in% method | "clr" %in% method) {
    print("CLR normalization")
    s <- Seurat::CreateAssayObject(protein_mat,assay="Protein")
    s <- Seurat::CreateSeuratObject(s,assay="Protein")
    s <- Seurat::NormalizeData(s, normalization.method = "CLR")
    SummarizedExperiment::assay(protein_sce, "CLR_norm") <- s@assays$Protein@data
  }
  if("SVD"%in%method|"svd"%in%method){
    
    temp_mat <- as.matrix(log10(t(protein_mat)+1))
    
    # calculates aspect Ratio
    aspect_ratio_beta_val <- (dim(temp_mat)[2])/(dim(temp_mat)[1])
    
    # Run SVD
    svd_protein <-svd(temp_mat)
    
    # Pull out the diagonal singularity terms
    truncated_singularity_values <- svd_protein$d

    # truncates the singularity terms based on the median marcenko-pastur formula 
    if(is.null(num_components_to_keep)){
      truncated_singularity_values[truncated_singularity_values<havok::optimal_SVHT_coef(aspect_ratio_beta_val,sigma_known = FALSE)*quantile(truncated_singularity_values)[3]]<-0
    }
    else if(num_components_to_keep<length(truncated_singularity_values)){
      truncated_singularity_values[num_components_to_keep:length(truncated_singularity_values)]=0
    }
    
    print(plot(truncated_singularity_values))
    # Reconstruction of the protein_matrix.
    protein_reconstructed <- svd_protein$u%*%diag(truncated_singularity_values)%*%t(svd_protein$v)
    rownames(protein_reconstructed)<- colnames(protein_mat)
    colnames(protein_reconstructed)<- rownames(protein_mat)
    SummarizedExperiment::assay(protein_sce, "SVD_norm")<-t(protein_reconstructed)
    
  }
  if(!all(method%in%c("dsb","CLR","SVD"))){
    print("method not available")
    break
  }
  
  protein_sce@metadata<-metadata
  protein_sce@colData<-S4Vectors::DataFrame(metadata%>%
                                            dplyr::filter(Cell%in%colnames(protein_mat)))
  SingleCellExperiment::altExp(sce, "Protein") <- protein_sce
  
  return(sce)
}
