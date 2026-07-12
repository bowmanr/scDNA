#' Normalize Protein Data
#' @description
#' This wrapper function provides available strategies for single cell protein to normalize the counts.
#'
#' @param sce the SingleCellExperiment object after running tapestri_h5_to_sce()
#' @param metadata this is droplet_metadata from extract_droplet_size() function.
#' @param method there are 3 methods: 1) a CLR approach, 2) a dsb approach, and 3) an experimental truncated-SVD approach using the marcenko-pastur formula
#' @param detect_IgG a flag to account for IgG controls in the panel.
#' @param num_components_to_keep Only used for the SVD approach on how many latent features to keep.
#' @param background_droplets a list of of names for the background empty droplets. Only used for dsb to disambiguate signal from noisy samples with background
#'
#' @return it returns the sce object with an additional alternative experiment
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @examples \dontrun{droplet_metadata<- extract_droplet_size(sce)
#' background_droplets<-droplet_metadata%>%
#' dplyr::filter(Droplet_type=="Empty")%>%
#' dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
#' dplyr::pull(Cell)
#' sce_subset<-normalize_protein_data(sce=sce,
#' metadata=droplet_metadata,
#' method=c("dsb"),
#' detect_IgG=TRUE,
#' background_droplets=background_droplets)
#' }
#'
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
                                    dplyr::select(tidyselect::all_of(background_droplets))

    rownames(empty_drops_matrix_input)<- rownames(protein_mat)

    if(detect_IgG){
        isotype <- grep("IgG",colnames(protein_mat),value=TRUE)
    }else{
        isotype <- FALSE
    }

    adt_norm <- dsb::DSBNormalizeProtein(
      # remove ambient protein noise reflected in counts from empty droplets
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
      truncated_singularity_values[truncated_singularity_values<havok::optimal_SVHT_coef(aspect_ratio_beta_val,sigma_known = FALSE)*stats::quantile(truncated_singularity_values)[3]]<-0
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
    stop("Method not available. Use one or more of: 'dsb', 'CLR', or 'SVD'.", call. = FALSE)
  }

  protein_sce@metadata<-metadata
  protein_sce@colData<-S4Vectors::DataFrame(metadata%>%
                                            dplyr::filter(Cell%in%colnames(protein_mat)))
  SingleCellExperiment::altExp(sce, "Protein") <- protein_sce

  return(sce)
}
