
#' Visualize Full Network 
#'
#' @param save_filename file location to save file. Currently, must end it with .html extension.
#' @param sce  SingleCellExperiment object containing NGT matrix for clone identification.
#' @return interactive for the Reinforcement Learning Trajectory
#' @export
#'
#' @examples
visualize_full_network<-function(sce,save_filename=NULL){
  
  RL<-sce@metadata$RL_info
  temp_traj<-RL%>%as.data.frame%>%
    dplyr::select(current_state,next_state,reward,action_type)%>%
    dplyr::rename(mutation_taken="action_type")
  new_net<-temp_traj

  new_net$current_state<- as.character(as.numeric(gsub("_", "",new_net$current_state)))
  new_net$next_state<- as.character(as.numeric(gsub("_", "",new_net$next_state)))
  #colnames(new_net)<-c("from","to","value","type")
  new_net<-new_net%>%dplyr::rename(from="current_state",
                            to="next_state",
                            value="reward",
                            type="mutation_taken")

  data_frame_clarity <- dplyr::inner_join(sce@metadata$Architecture, 
                                          as.data.frame(sce@metadata$Clones) %>% 
                                            dplyr::mutate(Count = ifelse(Group == "Other", n_Other, n_Complete)) %>% 
                                            dplyr::filter(Group =="Complete") %>% 
                                            dplyr::select(Clone, Count) %>% 
                                            distinct, by = "Clone") %>% dplyr::group_by(Clone)
  weights_unique_list <- distinct(data.frame(value = data_frame_clarity$Count/sum(data_frame_clarity$Count), 
                                             label = as.character(as.numeric(gsub("_", "", data_frame_clarity$Clone)))))
  
  
  
  
  data2<-visNetwork::toVisNetworkData(igraph::graph_from_data_frame(new_net))
  print(data2$nodes)
  data2$nodes<-data2$nodes%>%dplyr::inner_join(RL%>%as.data.frame%>%
                                          dplyr::select(next_state,observed_states)%>%
                                            distinct%>%
                                            dplyr::rename(label="next_state")%>%
                                            dplyr::rename(group="observed_states"),
                                        by="label")
  data2$nodes<-data2$nodes%>%dplyr::mutate(group=ifelse(group==1,
                                                 "observed",
                                                 "unobserved"))
  data2$nodes<-data2$nodes%>%dplyr::mutate(color=ifelse(group=="observed",
                                                 "darkred",
                                                 "grey"),
                                    font.size=10,
                                    scaling.min=4)
  print(data2$nodes)
  
  data2$nodes <-data2$nodes%>%
    dplyr::left_join(weights_unique_list,by="label")%>%
    dplyr::mutate(value=ifelse(is.na(value),0,value))
  data2$edges<-data2$edges %>%
    dplyr::mutate(color=dplyr::case_when(type=="none"~"blue",
                           type=="ADO"~"green",
                           type=="forward_ADO"~"magenta",
                           TRUE~"black"),
           dashes=dplyr::case_when(type=="none"~FALSE,
                            type=="ADO"~TRUE,
                            type=="forward_ADO"~TRUE,
                            TRUE~FALSE),
           label=type,
           title=type,
           arrows="to",
           arrowStrikethrough=F,
           scaling.label=F,
           font.size=dplyr::case_when(type=="ADO"~0,
                               type=="forward_ADO"~0,
                               TRUE~16),
           font.align="top",
           font.vadjust=(-30))
  
  
  if(is.null(save_filename)){
    visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
      visNetwork::visEdges(arrows = "to")%>%
      visNetwork::visLegend(useGroups = F,
                zoom = F,
                addNodes = data.frame(label=c("observed","unobserved"),
                                      shape=c("circle","circle"),
                                      color=c("darkred","grey")),
                addEdges = data.frame(label=c("ADO",
                                              "forward_ADO",
                                              "mutation",
                                              "no action"),
                                      color=c("green","magenta","black","blue"),
                                      dashes=c(TRUE,TRUE,FALSE,FALSE),
                                      font.align="top",
                                      font.size=20,
                                      font.vadjust=(-10)))%>%
      #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
      visNetwork::visLayout(randomSeed = 44)%>%  
      visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")%>%
      visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)#%>%visSave(file = "WT_to_match_clonograph.html",)
    
    
  }else{
    visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
      visNetwork::visEdges(arrows = "to")%>%
      visNetwork::visLegend(useGroups = F,
                zoom = F,
                addNodes = data.frame(label=c("observed","unobserved"),
                                      shape=c("circle","circle"),
                                      color=c("darkred","grey")),
                addEdges = data.frame(label=c("ADO",
                                              "forward_ADO",
                                              "mutation",
                                              "no action"),
                                      color=c("green","magenta","black","blue"),
                                      dashes=c(TRUE,TRUE,FALSE,FALSE),
                                      font.align="top",
                                      font.size=20,
                                      font.vadjust=(-10)))%>%
      #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
      visNetwork::visLayout(randomSeed = 44)%>%  
      visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")%>%
      visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)%>%
      visNetwork::visSave(file = save_filename)
  }
                        
}


#' BSCITE Style single trajectory Plot
#'
#' @param sce SingleCellExperiment object containing NGT matrix for clone identification.
#' @param trajectory Put in entire Trajectory you are interested in such as 'sce(at)metadata$Trajectories[[1]]'
#' @param save_filename file location to save file. Currently, must end it with .html extension.
#' @return
#' @export
#'
#' @examples
trajectory_of_interest_BSCITE_format<-function(sce,trajectory=sce@metadata$Trajectories[[1]],save_filename=NULL){
  temp_traj<-trajectory
  if(class(temp_traj)=="list"){
    new_net<-as.data.frame(do.call(rbind, temp_traj))%>%distinct
  }else{
    new_net<-temp_traj
  }
  
  
  net <- sce@metadata$RL_net
  RL <- sce@metadata$RL_info
  net <- igraph::graph_from_data_frame(data.frame(from = as.numeric(levels(RL$current_state))[RL$current_state], 
                                                  to = as.numeric(levels(RL$next_state))[RL$next_state], 
                                                  weight = abs(1/(RL$Q_values_normalized))), directed = TRUE)
  data_frame_clarity <- dplyr::inner_join(sce@metadata$Architecture, 
                                          as.data.frame(sce@metadata$Clones) %>% 
                                            dplyr::mutate(Count = ifelse(Group =="Other", n_Other, n_Complete)) %>% 
                                            dplyr::filter(Group == "Complete") %>% dplyr::select(Clone, Count) %>% 
                                            distinct, by = "Clone") %>% dplyr::group_by(Clone)
  weights_unique_list <- distinct(data.frame(value = data_frame_clarity$Count/sum(data_frame_clarity$Count), 
                                                    label = as.character(as.numeric(gsub("_", "", data_frame_clarity$Clone)))))
  
  
  countdf<- new_net %>%
    group_by(mutation_taken)%>%
    summarise(count=n())
  
  bscite_figure_df<-new_net %>% group_by(mutation_taken)%>%
    dplyr::left_join(countdf)%>%
    dplyr::mutate(r=row_number(),
           node_name = ifelse((mutation_taken!="forward_ADO"&mutation_taken!="ADO"),ifelse(r>1,paste0(mutation_taken,"_Hom"),paste0(mutation_taken,"_Het")),mutation_taken))
  
  
  bscite_figure_df$current_state<- as.character(as.numeric(gsub("_", "",bscite_figure_df$current_state)))
  bscite_figure_df$next_state<- as.character(as.numeric(gsub("_", "",bscite_figure_df$next_state)))
  
  
  
  bscite_figure_df<-bscite_figure_df%>%dplyr::rename(from="current_state",
                                              to="next_state",
                                              value="reward",
                                              type="node_name")
  
  data_bscite<-toVisNetworkData(igraph::graph_from_data_frame(bscite_figure_df))
  
  data_bscite$nodes<-data_bscite$nodes%>%
    dplyr::inner_join(RL%>%                              
                 dplyr::select(next_state,observed_states)%>%
                 distinct%>%                                      
                 dplyr::rename(label="next_state")%>%
                 dplyr::rename(group="observed_states"),
               by="label")
  data_bscite$nodes <-data_bscite$nodes%>%
    dplyr::left_join(weights_unique_list,by="label")%>%
    dplyr::mutate(value=ifelse(is.na(value),0,value))
  
  data_bscite$nodes<-data_bscite$nodes%>%
    dplyr::inner_join(data_bscite$edges%>%
                 dplyr::select(to,type)%>%
                 dplyr::rename(to="id"),by="id")%>%
    dplyr::select(-label)%>%
    dplyr::rename(type="label")
  
  data_bscite$nodes<-data_bscite$nodes%>%
    dplyr::mutate(font.vadjust=(0),
           shape ="ellipse",
           color.background="white",
           color.border = "black")#case_when(label== ~"black",
  #          label== ~""))
  
  
  data_bscite$edges<-data_bscite$edges %>%
    dplyr::mutate(color=dplyr::case_when(type=="none"~"blue",
                           type=="ADO"~"green",
                           type=="forward_ADO"~"magenta",
                           TRUE~"black"),
           dashes=dplyr::case_when(type=="none"~FALSE,
                            type=="ADO"~TRUE,
                            type=="forward_ADO"~TRUE,
                            TRUE~FALSE),
           #label=as.character(round(value,digits=4)),
           #title=as.character(round(value,digits=4)),
           arrows="to",
           arrowStrikethrough=F,
           scaling.label=F,
           font.size=dplyr::case_when(type=="ADO"~0,
                               type=="forward_ADO"~0,
                               TRUE~16),
           font.hadjust=(-50))%>%
    #font.vadjust=(-30),
    #font.align="top")%>%
    dplyr::rename(value="weight")
    
  if(is.null(save_filename)){
    visNetwork::visNetwork(nodes = data_bscite$nodes, edges = data_bscite$edges)%>%#,width = "100%",height = "100%")%>% 
      visNetwork::visEdges(arrows = "to")%>%
      #visLayout(randomSeed = 44)%>%  
      visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")%>%
      visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)
    
    }else{
      visNetwork::visNetwork(nodes = data_bscite$nodes, edges = data_bscite$edges)%>%#,width = "100%",height = "100%")%>% 
        visNetwork::visEdges(arrows = "to")%>%
        #visLayout(randomSeed = 44)%>%  
        visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")%>%
        visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)%>%
      visNetwork::visSave(file = save_filename)
    
    
  } 
}

#' Individual trajectory of interest Plots
#'
#' @param sce SingleCellExperiment object containing NGT matrix for clone identification.
#' @param trajectory Put in entire Trajectory you are interested in such as 'sce at metadata$Trajectories[[1]]'
#' @param save_filename file location to save file. Currently, must end it with .html extension.
#' @return
#' @export
#'
#' @examples
trajectory_of_interest_figure<-function(sce,trajectory=sce@metadata$Trajectories[[1]],save_filename=NULL){
  temp_traj<-trajectory
  tree_flag=FALSE
  if(class(temp_traj)=="list"){
    new_net<-as.data.frame(do.call(rbind, temp_traj))%>%distinct
  }else{
    new_net<-temp_traj
  }
  net <- sce@metadata$RL_net
  RL <- sce@metadata$RL_info
  net <- igraph::graph_from_data_frame(data.frame(from = as.numeric(levels(RL$current_state))[RL$current_state], 
                                                  to = as.numeric(levels(RL$next_state))[RL$next_state], 
                                                  weight = abs(1/(RL$Q_values_normalized))), directed = TRUE)
  data_frame_clarity <- dplyr::inner_join(sce@metadata$Architecture, 
                                          as.data.frame(sce@metadata$Clones) %>% 
                                            dplyr::mutate(Count = ifelse(Group =="Other", n_Other, n_Complete)) %>% 
                                            dplyr::filter(Group == "Complete") %>% dplyr::select(Clone, Count) %>% 
                                            distinct, by = "Clone") %>% dplyr::group_by(Clone)
  weights_unique_list <- distinct(data.frame(value = data_frame_clarity$Count/sum(data_frame_clarity$Count), 
                                             label = as.character(as.numeric(gsub("_", "", data_frame_clarity$Clone)))))
  
  new_net$current_state<- as.character(as.numeric(gsub("_", "",new_net$current_state)))
  new_net$next_state<- as.character(as.numeric(gsub("_", "",new_net$next_state)))
  #colnames(new_net)<-c("from","to","value","type")
  new_net<-new_net%>%dplyr::rename(from="current_state",
                            to="next_state",
                            value="reward",
                            type="mutation_taken")
  data2<-visNetwork::toVisNetworkData(igraph::graph_from_data_frame(new_net))
  
  
  data2$nodes<-data2$nodes%>%dplyr::inner_join(RL%>%
                                          dplyr::select(next_state,observed_states)%>%
                                               distinct%>%
                                          dplyr::rename(label="next_state")%>%
                                               dplyr::rename(group="observed_states"),
                                        by="label")
  data2$nodes<-data2$nodes%>%dplyr::mutate(group=ifelse(group==1,
                                                 "observed",
                                                 "unobserved"))
  data2$nodes<-data2$nodes%>%dplyr::mutate(color=ifelse(group=="observed",
                                                 "darkred",
                                                 "grey"),
                                    font.size=10,
                                    scaling.min=2)
  data2$nodes <-data2$nodes%>%
    dplyr::left_join(weights_unique_list,by="label")%>%
    dplyr::mutate(value=ifelse(is.na(value),0,value))
  data2$edges<-data2$edges %>%
    dplyr::mutate(color=dplyr::case_when(type=="none"~"blue",
                           type=="ADO"~"green",
                           type=="forward_ADO"~"magenta",
                           TRUE~"black"),
           dashes=dplyr::case_when(type=="none"~FALSE,
                            type=="ADO"~TRUE,
                            type=="forward_ADO"~TRUE,
                            TRUE~FALSE),
           label=type,
           title=type,
           arrows="to",
           arrowStrikethrough=F,
           scaling.label=F,
           font.size=case_when(type=="ADO"~0,
                               type=="forward_ADO"~0,
                               TRUE~16),
           font.align="top",
           font.vadjust=(-30))
  
  if(is.null(save_filename)){
    if(tree_flag){
    visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
    visNetwork::visEdges(arrows = "to")%>%
    visNetwork::visLegend(useGroups = F,
              zoom = F,
              addNodes = data.frame(label=c("observed","unobserved"),
                                    shape=c("circle","circle"),
                                    color=c("darkred","grey")),
              addEdges = data.frame(label=c("ADO",
                                            "forward_ADO",
                                            "mutation",
                                            "no action"),
                                    color=c("green","magenta","black","blue"),
                                    dashes=c(TRUE,TRUE,FALSE,FALSE),
                                    font.align="top",
                                    font.size=20,
                                    font.vadjust=(-10)))%>%
    #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
    visNetwork::visLayout(randomSeed = 44)%>%  
    visIgraphLayout(layout = "layout_with_sugiyama")%>%
    visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)#%>%visSave(file = "WT_to_match_clonograph.html",)
    }
    else{
      visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
    visNetwork::visEdges(arrows = "to")%>%
    visNetwork::visLegend(useGroups = F,
              zoom = F,
              addNodes = data.frame(label=c("observed","unobserved"),
                                    shape=c("circle","circle"),
                                    color=c("darkred","grey")),
              addEdges = data.frame(label=c("ADO",
                                            "forward_ADO",
                                            "mutation",
                                            "no action"),
                                    color=c("green","magenta","black","blue"),
                                    dashes=c(TRUE,TRUE,FALSE,FALSE),
                                    font.align="top",
                                    font.size=20,
                                    font.vadjust=(-10)))%>%
    #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
    visNetwork::visLayout(randomSeed = 44)%>%  
    #visIgraphLayout(layout = "layout_with_sugiyama")%>%
    visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)#%>%visSave(file = "WT_to_match_clonograph.html",)
    }
  }else{
    if(tree_flag){
    visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
      visNetwork::visEdges(arrows = "to")%>%
      visNetwork::visLegend(useGroups = F,
                            zoom = F,
                            addNodes = data.frame(label=c("observed","unobserved"),
                                                  shape=c("circle","circle"),
                                                  color=c("darkred","grey")),
                            addEdges = data.frame(label=c("ADO",
                                                          "forward_ADO",
                                                          "mutation",
                                                          "no action"),
                                                  color=c("green","magenta","black","blue"),
                                                  dashes=c(TRUE,TRUE,FALSE,FALSE),
                                                  font.align="top",
                                                  font.size=20,
                                                  font.vadjust=(-10)))%>%
      #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
      visNetwork::visLayout(randomSeed = 44)%>%  
      visIgraphLayout(layout = "layout_with_sugiyama")%>%
      visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)%>%
      visNetwork::visSave(file = save_filename)
    }
    else{
      visNetwork::visNetwork(nodes = data2$nodes, edges = data2$edges)%>%#,width = "100%",height = "100%")%>% 
      visNetwork::visEdges(arrows = "to")%>%
      visNetwork::visLegend(useGroups = F,
                            zoom = F,
                            addNodes = data.frame(label=c("observed","unobserved"),
                                                  shape=c("circle","circle"),
                                                  color=c("darkred","grey")),
                            addEdges = data.frame(label=c("ADO",
                                                          "forward_ADO",
                                                          "mutation",
                                                          "no action"),
                                                  color=c("green","magenta","black","blue"),
                                                  dashes=c(TRUE,TRUE,FALSE,FALSE),
                                                  font.align="top",
                                                  font.size=20,
                                                  font.vadjust=(-10)))%>%
      #visHierarchicalLayout(direction="UD",levelSeparation = 500)%>%
      visNetwork::visLayout(randomSeed = 44)%>%  
      #visIgraphLayout(layout = "layout_with_sugiyama")%>%
      visNetwork::visInteraction(dragNodes = TRUE, dragView = FALSE, zoomView = FALSE)%>%
      visNetwork::visSave(file = save_filename)
      }
    
  } 
}
      
