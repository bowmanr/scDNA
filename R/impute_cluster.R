#' Title
#'
#' @param sce 
#'
#' @return
#' @export
#'
#' @examples
impute_cluster<-function(sce,by="AF"){
  
  if(by=="AF"){
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  AF<-sce@assays@data$AF
  rownames(AF)<-rownames(sce)
  cluster_medians<-AF%>%
    tibble::rownames_to_column(var="Variant")%>%
    tidyr::pivot_longer(cols=-Variant,
                        values_to="AF",
                        names_to="Cell")%>%
    dplyr::full_join(SummarizedExperiment::colData(sce)%>%
                       data.frame%>%
                       dplyr::select(Cluster)%>%
                       tibble::rownames_to_column(var="Cell"),by="Cell")%>%
    dplyr::filter(!is.na(Cluster))%>%
    dplyr::group_by(Variant,Cluster)%>%
    dplyr::reframe(AF_med=median(AF))%>%
    dplyr::mutate(Cluster=factor(Cluster))%>%
    dplyr::arrange(Variant)%>%
    split(.,f=.$Cluster)
  
  missing_cell_AF<-sce[,is.na(SummarizedExperiment::colData(sce)$Cluster)]@assays@data$AF
  missing_cell_names<- colnames(missing_cell_AF)
  rownames(missing_cell_AF)<-rownames(sce)
  missing_cell_AF_long <-missing_cell_AF%>%
    tibble::rownames_to_column(var="Variant")%>%
    tidyr::pivot_longer(cols=-Variant,
                        values_to="AF",
                        names_to="Cell")%>%
    dplyr::group_by(Cell)
  
  cell_distances<-lapply(missing_cell_names,function(cellID){
    cell_AF_vec<-missing_cell_AF_long%>%dplyr::filter(Cell==cellID)%>%
      arrange(Variant) %>% pull(AF) %>%
      setNames(missing_cell_AF_long%>%dplyr::filter(Cell==cellID) %>% arrange(Variant) %>% pull(Variant))
    out<-lapply(seq_along(cluster_medians),function(clusterID){
      cluster_AF_vec<-cluster_medians[[clusterID]]%>%
        arrange(Variant) %>% pull(AF_med) %>%
        setNames(cluster_medians[[clusterID]] %>% arrange(Variant) %>% pull(Variant))
      return(euclidean(cell_AF_vec,cluster_AF_vec))
    })
    return(which(unlist(out)==min(unlist(out))))
  })
  new_metadata <- SummarizedExperiment::colData(sce)%>%
    data.frame%>%
    tibble::rownames_to_column(var="Cell")%>%
    dplyr::full_join( data.frame("Cell"=missing_cell_names,
                                 "imputed_cluster"=unlist(cell_distances)),
                      by="Cell")%>%
    dplyr::mutate(cluster_group=dplyr::case_when(
      is.na(imputed_cluster)~"Observed",
      TRUE~"Imputed"
    ))%>%
    dplyr::mutate(final_cluster=case_when(
      is.na(imputed_cluster)~Cluster,
      TRUE~imputed_cluster
    ))%>%
    dplyr::select(!imputed_cluster)
  return(add_cell_annotation(sce,new_metadata))
  } else if(by=="NGT"){
    hamming <- function(a, b) sum(a != b)
    NGT<-sce@assays@data$NGT
    rownames(NGT)<-rownames(sce)
    cluster_medians<-NGT%>%
      tibble::rownames_to_column(var="Variant")%>%
      tidyr::pivot_longer(cols=-Variant,
                          values_to="NGT",
                          names_to="Cell")%>%
      dplyr::full_join(SummarizedExperiment::colData(sce)%>%
                         data.frame%>%
                         dplyr::select(Cluster)%>%
                         tibble::rownames_to_column(var="Cell"),by="Cell")%>%
      dplyr::filter(!is.na(Cluster))%>%
      dplyr::group_by(Variant,Cluster)%>%
      dplyr::reframe(NGT_med=median(NGT))%>%
      dplyr::mutate(Cluster=factor(Cluster))%>%
      dplyr::arrange(Variant)%>%
      split(.,f=.$Cluster)
    
    missing_cell_NGT<-sce[,is.na(SummarizedExperiment::colData(sce)$Cluster)]@assays@data$NGT
    missing_cell_names<- colnames(missing_cell_NGT)
    rownames(missing_cell_NGT)<-rownames(sce)
    missing_cell_NGT_long <-missing_cell_NGT%>%
      tibble::rownames_to_column(var="Variant")%>%
      tidyr::pivot_longer(cols=-Variant,
                          values_to="NGT",
                          names_to="Cell")%>%
      dplyr::group_by(Cell)
    
    cell_distances<-lapply(missing_cell_names,function(cellID){
      cell_NGT_vec<-missing_cell_NGT_long%>%dplyr::filter(Cell==cellID)%>%
        arrange(Variant) %>% pull(NGT) %>%
        setNames(missing_cell_NGT_long%>%dplyr::filter(Cell==cellID) %>% arrange(Variant) %>% pull(Variant))
      out<-lapply(seq_along(cluster_medians),function(clusterID){
        cluster_NGT_vec<-cluster_medians[[clusterID]]%>%
          arrange(Variant) %>% pull(NGT_med) %>%
          setNames(cluster_medians[[clusterID]] %>% arrange(Variant) %>% pull(Variant))
        return(hamming(cell_NGT_vec,cluster_NGT_vec))
      })
      return(which(unlist(out)==min(unlist(out))))
    })
    new_metadata <- SummarizedExperiment::colData(sce)%>%
      data.frame%>%
      tibble::rownames_to_column(var="Cell")%>%
      dplyr::full_join( data.frame("Cell"=missing_cell_names,
                                   "imputed_cluster"=unlist(cell_distances)),
                        by="Cell")%>%
      dplyr::mutate(cluster_group=dplyr::case_when(
        is.na(imputed_cluster)~"Observed",
        TRUE~"Imputed"
      ))%>%
      dplyr::mutate(final_cluster=case_when(
        is.na(imputed_cluster)~Cluster,
        TRUE~imputed_cluster
      ))%>%
      dplyr::select(!imputed_cluster)
    return(add_cell_annotation(sce,new_metadata))
  } else {
    print("'by' not specified, provide AF or NGT")
  }
}