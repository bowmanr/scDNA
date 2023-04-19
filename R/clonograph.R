
#' Plotting clonographs
#'
#' @param final_sample_summary
#' @param complete_only logical to decide whether to plot only complete cells or also low quality cells, default=FALSE
#' @param color_pal brewer.pal pallete selection.
#' @importFrom magrittr %<>%
#' @return
#' @export
#'
#' @examples
clonograph<-function(final_sample_summary,
                     complete_only=FALSE,
                     color_pal="Reds"){

# Extract out the sample of interest
consolidated_clonal_abundance <-final_sample_summary$Clones%>%
                                    dplyr::group_by(Clone,Group)%>%
                                    dplyr::mutate(AF_med=mean(AF_med),
                                                     DP_med=mean(DP_med),
                                                     GQ_med=mean(GQ_med))%>%
                                    dplyr::select(!c(variants,Count))%>%
                                    dplyr::distinct()%>%
                                    dplyr::rowwise()%>%
                                    dplyr::mutate(Count=ifelse(Group=="Other",n_Other,n_Complete))%>%
                                    dplyr::ungroup()

if(complete_only==TRUE){
  consolidated_clonal_abundance <-consolidated_clonal_abundance%>%dplyr::filter(.data$Group=="Complete")
  }

consolidated_clonal_abundance%<>%dplyr::arrange(.data$Count)
clonal_architecture <-final_sample_summary$Architecture
mutant_order <-setdiff(colnames(final_sample_summary$NGT),c("Cell","Clone","Group"))

# Ensure the order of the clone abundance and clone architecture are the same.
clonal_architecture$Clone <- factor(clonal_architecture$Clone, levels=unique(rev(consolidated_clonal_abundance$Clone)))
consolidated_clonal_abundance$Clone <- factor(consolidated_clonal_abundance$Clone, levels=levels(clonal_architecture$Clone))
consolidated_clonal_abundance$Group <- factor(consolidated_clonal_abundance$Group, levels=c("Complete","Other"))

# Generate clonal abundance barplot
gg_clonal_barplot <- ggplot(data=consolidated_clonal_abundance, aes(x=Clone, y=Count,fill=Group)) +
  geom_bar(fun = "identity", stat = "summary",position = "stack")+
  theme_classic(base_size=7)+
  scale_y_continuous(expand=c(0.01,0))+
  #ylim() +
  ylab("Cell Count")+
 {if(complete_only==FALSE) {geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2)}}+
  scale_fill_manual(values=c("Other"="Grey70",
                              "Complete"=RColorBrewer::brewer.pal(5,tidyselect::all_of(color_pal))[5])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x =element_blank(),
        legend.position = "right",
        plot.margin=unit(c(0,0,0,0),"cm"))

# Generate mutation heatmap
gg_heatmap <- ggplot(data=clonal_architecture,
                     aes(x=Clone,y=AA,fill=Genotype))+
  geom_tile() +
  scale_fill_manual(values=c("WT"=brewer.pal(7,tidyselect::all_of(color_pal))[1],
                             "Heterozygous"=brewer.pal(7,tidyselect::all_of(color_pal))[3],
                             "Homozygous"=brewer.pal(7,tidyselect::all_of(color_pal))[6],
                             "Unknown"="grey50"),name="Genotype")+
  theme_classic(base_size=7) +
  ylab("Mutation")+
  scale_y_discrete(limits = rev(levels(clonal_architecture$AA)))+
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))

gg_QC_heatmap <- ggplot(data=consolidated_clonal_abundance,
                     aes(x=Clone, y=Group, fill=ADO_med))+
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0.1,
                                                    rev=FALSE,na.value = "grey80",limits=c(0,0.25))+
  theme_classic(base_size=7) +
  theme(legend.position = "right", legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))


gg_QC_heatmap_GQ <- ggplot(data=consolidated_clonal_abundance,
                        aes(x=Clone, y=Group, fill=GQ_med))+
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 30,
                                               rev=TRUE,na.value = "grey80",limits=c(0,100))+
  theme_classic(base_size=7) +
  theme(legend.position = "right", legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))

gg_QC_heatmap_DP <- ggplot(data=consolidated_clonal_abundance,
                           aes(x=Clone, y=Group, fill=DP_med))+
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 10,
                                               rev=TRUE,na.value = "grey80",limits=c(0,40))+
  theme_classic(base_size=7) +
  theme(legend.position = "right", legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))

# Put it all together
return(plot_grid(gg_clonal_barplot,gg_heatmap,gg_QC_heatmap,gg_QC_heatmap_GQ,gg_QC_heatmap_DP,ncol=1,align="v",axis="lr",rel_heights = c(1,0.5,0.25,0.25,0.25)))
}
