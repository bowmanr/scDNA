#' Validate SCE
#' @description
#' This function returns a report of what analyses have been conducted and still are available to be completed.
#'
#' @param sce The current singleCellExperiment object.
#'
#' @return it outputs a status report directly to the console about processes that have been completed or
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'     validate_sce(sce)
#' }
validate_sce<-function(sce){


  # enumerate clone check
  if(c("Architecture")%in%names(sce@metadata)){
    message("Clonal Architecture found: enumerate_clones() completed")
  }else{
    warning("Clonal Architecture NOT found: please proceed with enumerate_clones()")
  }

  # compute clone statistics check
  if(c("sample_stats")%in%names(sce@metadata)){
    message("Sample Stats found: compute_clone_statistics() has been completed")
  }else{
    warning("Sample Stats NOT found: please proceed with compute_clone_statistics()")
  }

  # CNV read from tapestri_h5 check
  if(c("CNV")%in%SingleCellExperiment::altExpNames(sce)){
    message("CNV data found: tapestri_h5_to_sce() completed")
    # ploidy from readDNA_CN_H5 check via compute_clone_statistics
    if(c("ploidy")%in%names(sce@metadata)){
      message("ploidy data has been found: compute_clone_statistics() and readDNA_CN_H5() have been completed")
    }else{
      warning("ploidy NOT found: please run either compute_clone_statistics(skip_ploidy=FALSE), or readDNA_CN_H5(). Note: enumerate_clones() needs to be ran prior compute_clone_statistics() function but readDNA_CN_H5() does not.")
    }
  }
  else{
    warning("CNV data NOT found: please run tapestri_h5_to_sce()")
  }

  # Protein checks
  if(c("Protein")%in%SingleCellExperiment::altExpNames(sce)){

    # now check if we have run normalize_protein_data
    if(any(c("DSB_norm","CLR_norm","SVD_norm")%in%SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce,"Protein")))){
      message("Normalized protein assays found: extract_droplet_size() and normalize_protein_data() have been completed")

    }else{
      warning("Normalized protein assays NOT found: please proceed with extract_droplet_size() and normalize_protein_data()")
    }
  }else{
      #check if it exists. If it exists, then we need to recommend all 3
      skip <- tryCatch( rhdf5::h5read(file = sce@metadata$file,
                                      name = "/assays/protein_read_counts/layers/read_counts")%>%
                          nrow()%>%{.>0},
                        error = function(e) {
                          message("'Protein' dataset not found in h5. No protein analysis needed")
                          return(FALSE)
                          })
      if(skip){
        warning("Protein dataset in h5 present. but Protein Data in SCE NOT Found: Please run tapestri_h5_to_sce(protein=T),extract_droplet_size(), and normalize_protein_data() ")
      }
  }



  if(c("cell_label_confidence")%in%names(SummarizedExperiment::colData(sce))){
    message("Cell Confidence found: cell_confidence_labeling() completed.")

  }else{
    skip <- tryCatch( rhdf5::h5read(file = sce@metadata$file,
                                    name = "/assays/protein_read_counts/layers/read_counts")%>%
                        nrow()%>%{.>0},
                      error = function(e) {
                        message("'Protein' dataset not found in h5. No protein analysis needed")
                        return(FALSE)
                      })
    if(skip){
      warning("Protein dataset in h5 present, Cell Confidence NOT Found: Please run cell_confidence_labeling()")
    }

  }

  # RL check
  if(c("RL_info")%in%names(sce@metadata)){
    message("RL_info found: trajectory_analysis() completed")
  }else{
    warning("RL_info and trajectories NOT found: please proceed with trajectory_analysis()")
  }



}
