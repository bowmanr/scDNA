clonograph_new<-function (sce, complete_only = FALSE, color_pal = "Reds", QC_stats = FALSE) 
{
    
    consolidated_clonal_abundance <- sce@metadata$Clonal_Abundance %>% 
      dplyr::arrange(.data$Count)%>%
      tidyr::pivot_longer(cols=c(n_Other,n_Complete),names_to = "Group",values_to = "Total_Count")%>%
      dplyr::mutate(Group=ifelse(Group=="n_Other","Other","Complete"))
  if (complete_only == TRUE) {
    consolidated_clonal_abundance <- consolidated_clonal_abundance %>% 
      dplyr::filter(Group=="Complete")%>%
      dplyr::mutate(Count=Total_Count)%>%
      dplyr::arrange(.data$Count) %>% 
      dplyr::select(-LCI, -UCI) %>% 
      dplyr::rename(LCI = Complete_LCI, 
                    UCI = Complete_UCI)%>%
      dplyr::filter(!is.na(LCI))
  } else {
    consolidated_clonal_abundance <- consolidated_clonal_abundance %>% 
    dplyr::mutate(Count=Total_Count)%>%
    dplyr::arrange(.data$Count)
  }
    
  
  clonal_architecture <- sce@metadata$Architecture
  mutant_order <- setdiff(colnames(sce@metadata$NGT), c("Cell", "Clone", "Group"))
  
  clonal_architecture$Clone <- factor(clonal_architecture$Clone, 
                                      levels = unique(rev(consolidated_clonal_abundance$Clone)))
  consolidated_clonal_abundance$Clone <- factor(consolidated_clonal_abundance$Clone, 
                                                levels = levels(clonal_architecture$Clone))
  consolidated_clonal_abundance$Group <- factor(consolidated_clonal_abundance$Group, 
                                                levels = c("Complete", "Other"))
  gg_clonal_barplot <- ggplot(data = consolidated_clonal_abundance, 
                              aes(x = Clone, y = Count, fill = Group)) + geom_bar(fun = "identity", 
                                                                                  stat = "summary", position = "stack") + theme_classic(base_size = 7) + 
    scale_y_continuous(expand = c(0.01, 0)) + ylab("Cell Count") + 
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) + 
    scale_fill_manual(values = c(Other = "Grey70", Complete = RColorBrewer::brewer.pal(n = 5, 
                                                                                       name = color_pal)[5])) + theme(axis.title.x = element_blank(), 
                                                                                                                      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                                                      axis.line.x = element_blank(), legend.position = "right", 
                                                                                                                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
  gg_heatmap <- ggplot(data = clonal_architecture, aes(x = Clone, 
                                                       y = final_annot, fill = Genotype)) + geom_tile() + scale_fill_manual(values = c(WT = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                     color_pal)[1], Heterozygous = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                                                                            color_pal)[3], Homozygous = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                                                                                                                                 color_pal)[6], Unknown = "grey50"), name = "Genotype") + 
    theme_classic(base_size = 7) + ylab("Mutation") + scale_y_discrete(limits = rev(mutant_order)) + 
    theme(legend.position = "right", legend.direction = "vertical", 
          axis.text.x = element_blank(), axis.line = element_blank(), 
          axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
    return(cowplot::plot_grid(gg_clonal_barplot, gg_heatmap, 
                              ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 
                                                                                  0.5)))
}


