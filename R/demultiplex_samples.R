#' Demultiplex wrapper function to handle splitting samples
#'
#' @param sce a SingleCellExperiment object containing the Clones,Architecture, and NGT matrix
#' @param sensitivity_threshold vector of sensitivity thresholds for rows (variants) and columns (cells)
#' @param expected_samples the number of clusters or samples you expect
#'
#' @return a dataframe of the cells assigned to a specific cluster.
#' @export
#'
#' @examples
demultiplex_samples<-function(sce,sensitivity_threshold=c(0.01,0.0001),expected_samples=5){
  
  AF_out<-scDNA::optimize_matrix(sce,sensitivity_threshold = c(0.01,0.0001))
  AF_complete<-t(sce[colnames(AF_out),rownames(AF_out)]@assays@data$AF)
  kmeans_AF <- kmeans(AF_complete, centers = expected_samples, nstart = 30) #centers set based on number of samples in the multiplex
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
  
  umap_df_final %>% ggplot(aes(x = UMAP1, y = UMAP2))+
    geom_point(aes( color = Cluster), alpha = .2, size =2, show.legend = TRUE)+
    scale_color_manual(values = rev(c(tol(5))))+
    theme_classic()+
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size =10),
          plot.title = element_text(hjust = 0.5, size = 14))+
    labs(x = "UMAP1",y = "UMAP2")+
    facet_grid(~sample)
  superheat(AF_complete, 
            heat.pal = kovesi.rainbow_bgyrm_35_85_c71(n=100),
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
  
  cell_clust_ids<-data.frame(cell_names=colData(sce)@rownames,
                             final_cluster=colData(sce)$final_cluster)
  
  return(cell_clust_ids)
}


