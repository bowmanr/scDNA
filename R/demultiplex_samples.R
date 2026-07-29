#' Demultiplex wrapper function to handle splitting samples
#' @description
#' This function takes a singleCellExperiment and applies an iterative optimization strategy to identify the best full rank matrix of variant allele frequency and cells. Once this is established it is a kmeans clustering is performed to identify the number of expected samples. Cells not used for the kmeans are then imputed for their cluster.
#'
#' @param sce a SingleCellExperiment object containing the Clones,Architecture, and NGT matrix
#' @param sensitivity_threshold vector of sensitivity thresholds for rows (variants) and columns (cells)
#' @param expected_samples the number of clusters or samples you expect. We recommend adding one additional group to account for doublets and noise.
#'
#' @return a dataframe of the cells assigned to a specific cluster.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' cell_clust_df<-demultiplex_samples(sce,expected_samples=5)
#' }
demultiplex_samples<-function(sce,sensitivity_threshold=c(0.01,0.0001),expected_samples=5){

  AF_out<-optimize_matrix(sce,sensitivity_threshold)
  AF_complete<-t(sce[colnames(AF_out),rownames(AF_out)]@assays@data$AF)
  kmeans_AF <- stats::kmeans(AF_complete, centers = expected_samples, nstart = 30) 
  umap_SNPs <-umap::umap(AF_complete, n_neighbors = sqrt(nrow(AF_complete)))

  umap_df <- umap_SNPs$layout %>%
    data.frame() %>%
    dplyr::mutate(Cluster=as.character(kmeans_AF$cluster))%>%
    dplyr::select(UMAP1=X1,UMAP2=X2,Cluster)%>%
    tibble::rownames_to_column(var = "Cell")

  umap_df_final <- SummarizedExperiment::colData(sce)%>%
    data.frame()%>%
    tibble::rownames_to_column(var = "Cell")%>%
    dplyr::select(Cell,sample)%>%
    dplyr::inner_join(umap_df,by="Cell")

  umap_df_final %>% ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2))+
    ggplot2::geom_point(ggplot2::aes( color = Cluster), alpha = .2, size =2, show.legend = TRUE)+
    ggplot2::scale_color_manual(values = rev(pals::tol(5)))+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12),
          axis.title.y = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(size = 10),
          axis.text.y = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size =10),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14))+
    ggplot2::labs(x = "UMAP1",y = "UMAP2")+
    ggplot2::facet_grid(~sample)
  superheat::superheat(AF_complete,
            heat.pal = pals::kovesi.rainbow_bgyrm_35_85_c71(n=100),
            scale = FALSE,
            heat.pal.values = seq(from=0,to=1,by=0.01),
            bottom.label.text.angle = 90,
            membership.rows = umap_df$Cluster,
            bottom.label.size = 0.2,
            bottom.label.text.size = 0,
            left.label.text.size = 2.8,
            left.label.size = .1,
            pretty.order.rows = TRUE,
            pretty.order.cols = TRUE,
            smooth.heat = F)
  
  table(umap_df_final$sample,umap_df_final$Cluster)
  
  sce<-add_cell_annotation(sce,umap_df_final)
  sce<-impute_cluster(sce,by="AF")

  cell_clust_ids<-data.frame(cell_names=SummarizedExperiment::colData(sce)@rownames,
                             final_cluster=SummarizedExperiment::colData(sce)$final_cluster)

  return(cell_clust_ids)
}


