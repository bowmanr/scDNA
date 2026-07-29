#' Get Stain Index for cells
#' @description
#' This function obtains a StainIndex metric used in flow cytometry. It obtains the median protein and DNA total reads per cell for empty and Cell droplets. It then divides this by the empty cells standard deviation.
#'
#' @param droplet_metadata this is the droplet metadata after running extract_droplet_size()
#'
#' @return returns a droplet_metadata dataframe with the additional columns for stain index of proteins and DNA
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' droplet_metadata<- extract_droplet_size(sce)
#' droplet_metadata<-GetStainIndex(droplet_metadata)
#' }
GetStainIndex<-function(droplet_metadata){
  empty_med_vec<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Empty")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::mutate(dna_median=median(dna_size),
           pro_median=median(protein_size))%>%
    dplyr::select(dna_median,pro_median)%>%dplyr::distinct()
  cell_med_vec<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Cell")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::mutate(dna_median=median(dna_size),
           pro_median=median(protein_size))%>%
    dplyr::select(dna_median,pro_median)%>%dplyr::distinct()
  sd.mat3<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Empty")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::summarise(cov_mat=cov(.))%>%as.matrix
  cell_ind_vec<-droplet_metadata%>%
    dplyr::select(dna_size,protein_size)

  droplet_metadata<-cbind(droplet_metadata,
                          as.matrix(sweep(as.matrix(cell_ind_vec),2,as.matrix(empty_med_vec)))%*%solve(2*(diag(diag(sd.mat3))**0.5)))

  colnames(droplet_metadata)[length(colnames(droplet_metadata))-1]<-"Stain_index_DNA_size"
  colnames(droplet_metadata)[length(colnames(droplet_metadata))]<-"Stain_index_protein_size"

  return(droplet_metadata)

}
#' Get membership quality of cells vs empty
#' @description
#' This metric performs Expectation Maximization (EM) between two droplet classes, Empty and Cell. This measure is helfpul as a Quality Control check to ensure high quality droplets are called cells for downstream analysis.
#'
#' @param droplet_metadata this is the droplet metadata after running extract_droplet_size()
#'
#' @return returns a droplet_metadata dataframe with the additional columns for posterior classification of Cell vs Empty based on total DNA and Protein Content. It also contains a column for the final assignment call.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' droplet_metadata<- extract_droplet_size(sce)
#' droplet_metadata<-GetMembership(droplet_metadata)
#' }
#'
GetMembership<-function(droplet_metadata){
  Cell_cov_mat<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Cell")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::summarise(cov_mat=cov(.))%>%as.matrix
  cell_mu_vec<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Cell")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::mutate(dna_mean=mean(dna_size),
           pro_mean=mean(protein_size))%>%
    dplyr::select(dna_mean,pro_mean)%>%dplyr::distinct()

  Empty_cov_mat<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Empty")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::summarise(cov_mat=cov(.))%>%as.matrix
  empty_mu_vec<-droplet_metadata%>%
    dplyr::filter(Droplet_type=="Empty")%>%
    dplyr::select(dna_size,protein_size)%>%
    dplyr::mutate(dna_mean=mean(dna_size),
           pro_mean=mean(protein_size))%>%
    dplyr::select(dna_mean,pro_mean)%>%dplyr::distinct()

  mu.vector1 <-t(cell_mu_vec%>%as.matrix)
  mu.vector1 <-mu.vector1[1:2]
  mu.vector2<-t(empty_mu_vec%>%as.matrix)
  mu.vector2 <-mu.vector2[1:2]

  sd.mat1 <- Cell_cov_mat
  sd.mat2 <-Empty_cov_mat

  alpha.vector1 <-1
  alpha.vector2 <-1
  comp1.prod=matrix(data=NA_real_,nrow = dim(droplet_metadata)[1],ncol=1)
  comp2.prod=matrix(data=NA_real_,nrow = dim(droplet_metadata)[1],ncol=1)
  for(iter in 1:dim(droplet_metadata)[1]){
    temp_df <-droplet_metadata[iter,c(2,5)]
    comp1.prod[iter]<-mvtnorm::dmvnorm(x=temp_df,
                                       mean=mu.vector1,
                                       sigma=(sd.mat1))*alpha.vector1[1]

    comp2.prod[iter] <- mvtnorm::dmvnorm(x=temp_df,
                                         mean=mu.vector2,
                                         sigma=(sd.mat2))*alpha.vector2[1]

  }
  droplet_metadata<-cbind(droplet_metadata,comp1.prod,comp2.prod)
  droplet_metadata<-droplet_metadata%>%
    dplyr::mutate(sum.of.comps=comp1.prod+comp2.prod,
           comp1.post= comp1.prod/sum.of.comps,
           comp2.post= comp2.prod/sum.of.comps)

  droplet_metadata<-droplet_metadata%>%
    dplyr::mutate(positive_cell_calls=ifelse(comp1.post>comp2.post,"Cell","Empty"))
  return(droplet_metadata)
}
