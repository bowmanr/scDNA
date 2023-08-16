mutant_clone_fractional_plot<-function(x_summary){
  
  
  clones<-as.data.frame(x_summary$Clones)%>%
    select(Clone)
  # Compute the relative clone size
  Clone_size <-  as.data.frame(x_summary$Clones)%>%
    transmute(Clone_size=Count/sum(Count/3))%>%    
    pull(Clone_size)
  names(Clone_size) <-t(clones[,"Clone"])
                         
  mutations <-colnames(as.data.frame(x_summary$NGT%>%
                    select(!Clone &!Group &!Cell)))
                         
  dominant_list_clone <- as.data.frame(x_summary$Clones)%>%
   arrange(desc(Count))%>%
   select(Clone)
 dominant_clone<-dominant_list_clone$Clone[1] 
  
  VAFs <-as.data.frame(x_summary$NGT)%>%
  select(!Clone &!Group &!Cell)%>%
  summarise_all(funs(mean(.)/2))
  
  
  mut_mat <-as.data.frame(x_summary$Clones)%>%
  select(Clone)%>%
  separate(col=Clone,into=mutations,sep="_")
  # creates a composite data frame and turn it into long form
  clone_size_by_gene<-data.frame(clones,Clone_size,mut_mat,check.names=FALSE)%>%
  pivot_longer(cols=all_of(mutations),names_to="Variant", values_to="Genotype")%>%
  filter(Genotype!=0)%>% # remove WT entries
  separate(col=Variant, into="Gene",extra="drop",sep="\\.|_",remove=FALSE)%>% # For later useage in plotting
  group_by(Variant)%>% 
  filter(Clone_size==max(Clone_size))%>% #identify largest clone
  mutate(Clonality=case_when(
   Clone==`dominant_clone`~"Dominant",
   Clone!=`dominant_clone`~"Subclone"))%>% #label clones
  inner_join(data.frame("VAF"=t(VAFs),
                       "Variant"=names(VAFs))) # merge with bulk VAF info
  
  
  tally_set <-data.frame(table(clone_size_by_gene$Gene,
                            clone_size_by_gene$Clonality))%>%
  pivot_wider(names_from=Var2,values_from=Freq)%>%
  mutate(Ratio=Dominant/(Subclone+Dominant))%>% #calculate the dominant ratio
  arrange(Ratio) 
  
  # the following are for the plots
  # For plotting purposes establish order of the y axis
  clone_size_by_gene$Gene <- factor(clone_size_by_gene$Gene, levels=tally_set$Var1)
  
  # Linde and I spent too much time picking the exact shade of red we wanted for this paper....
  color_red<-brewer.pal(5,"Reds")[5]
  
  # For plotting purposes establish order of stacked bars
  clone_size_by_gene$Clonality<-factor(clone_size_by_gene$Clonality,levels=c("Subclone","Dominant"))
  
  #plot the data
  ggA<-ggplot(tally(clone_size_by_gene%>%group_by(Gene,Clonality)),
           aes(x=factor(Gene),fill=Clonality,y=n,label=n))+
  guides(fill=FALSE,color=FALSE)+
  scale_y_continuous( expand = c(0, 0.0))+ #removes white space near the axis of the bars
  geom_bar(stat="identity",position="fill")+
  xlab("")+coord_flip()+
  scale_fill_manual(values=c("Dominant"=color_red,
                            "Subclone"="grey80"))+
  ylab("Fraction of mutant clones \n with mutation in dominant clone")+
  theme_bw(base_size=8)+theme(legend.position = "bottom")
  
  ## Need to put this somewhere up top
  
  ggB<-ggplot(clone_size_by_gene, 
           aes(y=Clone_size, x=Gene, fill=Gene)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(aes(color=Clonality,group=Clonality), 
            position = position_jitterdodge(), size=0.3)+
  scale_fill_manual(values=tol.rainbow(n=length(levels(clone_size_by_gene$Gene))))+
  scale_color_manual(values=c("Dominant"=color_red,
                             "Subclone"="grey20"))+
  coord_flip()+
  theme_bw(base_size=8)+guides(fill=FALSE,color=FALSE)+
  theme(axis.text.y = element_blank(),
       axis.ticks.y  = element_blank(),
       axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0.05)) +
  ylab("Fraction of cells \n in largest mutant clone")+
  theme(legend.position = "bottom")
  
  spacer <- plot_grid(NULL) # plot looks better with a little spacer
  
  plot_grid(ggA,spacer,ggB,align="h",axis="tb",
         ncol=3,rel_widths=c(1,0.05,1))

  
}