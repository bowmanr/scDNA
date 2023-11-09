# some things that need to be fixed a bit more I think!
# 1) I think people will flip start and goal state, so we will need to check which 
# one is the lower state and swap them (since we can only move forward). 
# 2) how should people input this? 000 vs 1000 vs 221 will not be intuitive for 
# the biologists I think. Is there a way you prefer where they give us strings 
# instead then we convert it to our 100 representation?

#' Title
#'
#' @param sce 
#' @importFrom igraph V
#' @return
#' @export
#'
#' @examples
visualize_full_network<-function(sce){
  net<-sce@metadata$RL_net
  RL <- sce@metadata$RL_info
    visNet<-visIgraph(net)%>%
    #visNetwork(nodes=visNet$nodes,edges=visNet$edges,width="100%")%>%
    visNodes(color=ifelse(names(V(net))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey"),
             #label=paste0("s",seq(0:(length(V(net))-1)))
    )%>%
    visEdges(color="black")%>%
    visIgraphLayout(layout="layout_with_sugiyama")
  visNet$x$options$nodes$font=list(size=48)
  visNet$x$nodes$color<-ifelse(names(V(net))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  visNet$x$nodes$label <-paste0("s",seq(0:(length(V(net))-1)))
  visNet$x$options$width="100%"
  visNet$x$nodes$title<-c("TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WW 0",
                          "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WM 1",
                          "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: MM 2",
                          
                          "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WW 0",
                          "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WM 1",
                          "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: MM 2",
                          
                          "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WW 0",
                          "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WM 1",
                          "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: MM 2",
                          
                          "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WW 0",
                          "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WM 1",
                          "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: MM 2",
                          
                          "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WW 0",
                          "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WM 1",
                          "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: MM 2",
                          
                          "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WW 0",
                          "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WM 1",
                          "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: MM 2",
                          
                          "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WW 0",
                          "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WM 1",
                          "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: MM 2",
                          
                          "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WW 0",
                          "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WM 1",
                          "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: MM 2",
                          
                          "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WW 0",
                          "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WM 1",
                          "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: MM 2")
  return(visNet)
}


#' Title
#'
#' @param sce 
#'
#' @return
#' @export
#'
#' @examples
visualize_WT_dominant_clone<-function(sce){
  net<-sce@metadata$RL_net
  RL <- sce@metadata$RL_info
  Total_hover_list <- c("TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: MM 2")
   Total_label_list<-paste0("s",seq(0:(length(V(net))-1)))
  
  WT_to_dom_clone <-igraph::shortest_paths(net,from="0",to=RL$next_state[which(RL$reward==max(RL$reward))[1]][1],algorithm ="dijkstra")$vpath[[1]]
  V(net)[WT_to_dom_clone]$color <-6
  g <- induced.subgraph(net,V(net)[WT_to_dom_clone])
  
  visNet<-visIgraph(g)%>%
    #visNetwork(nodes=visNet$nodes,edges=visNet$edges,width="100%")%>%
    visNodes(color=ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey"),
             #label=paste0("s",seq(0:(length(V(net))-1)))
    )%>%
    visEdges(color="black")%>%
    visIgraphLayout(layout="layout_with_sugiyama")
  visNet$x$options$nodes$font=list(size=48)
  visNet$x$options$width="100%"
  visNet$x$nodes$color<-ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  
  visNet$x$nodes$label <-Total_label_list[WT_to_dom_clone]
  visNet$x$nodes$title<-Total_hover_list[WT_to_dom_clone]
  return(visNet)
  
  
}

#' Title
#'
#' @param sce 
#'
#' @return
#' @export
#'
#' @examples
visualize_all_WT_dominant_clone<-function(sce){
  net<-sce@metadata$RL_net
  RL <- sce@metadata$RL_info
  Total_hover_list <- c("TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: MM 2")
  Total_label_list<-paste0("s",seq(0:(length(V(net))-1)))
  
  All_WT_to_dom<-all_simple_paths(net,from="0",to=RL$next_state[which(RL$reward==max(RL$reward))[1]][1])
  E(net)[unlist(All_WT_to_dom)]$weight = log(E(net)[unlist(All_WT_to_dom)]$weight)
  g <-induced.subgraph(net,E(net)[unlist(All_WT_to_dom)])
  
  visNet<-visIgraph(g)%>%
    #visNetwork(nodes=visNet$nodes,edges=visNet$edges,width="100%")%>%
    visNodes(color=ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey"),
             #label=paste0("s",seq(0:(length(V(net))-1)))
    )%>%
    visEdges(color="black")%>%
    visIgraphLayout(layout="layout_with_sugiyama")
  visNet$x$options$nodes$font=list(size=48)
  visNet$x$options$width="100%"
  visNet$x$nodes$color<-ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  
  visNet$x$nodes$label <-Total_label_list[unique(unlist(All_WT_to_dom))[order(unique(unlist(All_WT_to_dom)),decreasing = FALSE)]]
  visNet$x$nodes$title<-Total_hover_list[unique(unlist(All_WT_to_dom))[order(unique(unlist(All_WT_to_dom)),decreasing = FALSE)]]
  return(visNet)
  
  
}

#' Title
#'
#' @param sce 
#' @param start_name 
#' @param goal_name 
#'
#' @return
#' @export
#'
#' @examples
visualize_any_optimal_path <-function(sce,start_name,goal_name){
    net<-sce@metadata$RL_net
    RL <- sce@metadata$RL_info
    start_name<-as.character(as.numeric(gsub("_", "",start_name)))
    goal_name<-as.character(as.numeric(gsub("_", "",goal_name)))
    
  Total_hover_list <- c("TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: MM 2")
  Total_label_list<-paste0("s",seq(0:(length(V(net))-1)))
  clone_to_goal_clone <-igraph::shortest_paths(net,from=start_name,to=goal_name,algorithm = "dijkstra")$vpath[[1]]
  V(net)[clone_to_goal_clone]$color <-6
  g<- induced.subgraph(net,V(net)[clone_to_goal_clone])
  visNet<-visIgraph(g)%>%
    #visNetwork(nodes=visNet$nodes,edges=visNet$edges,width="100%")%>%
    visNodes(color=ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey"),
             #label=paste0("s",seq(0:(length(V(net))-1)))
    )%>%
    visEdges(color="black")%>%
    visIgraphLayout(layout="layout_with_sugiyama")
  visNet$x$options$nodes$font=list(size=48)
  visNet$x$options$width="100%"
  visNet$x$nodes$color<-ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  
  visNet$x$nodes$label <-Total_label_list[clone_to_goal_clone]
  visNet$x$nodes$title<-Total_hover_list[clone_to_goal_clone]
  return(visNet)
}

#' Title
#'
#' @param sce 
#' @param given_state 
#'
#' @return
#' @export
#'
#' @examples
match_clonal_graph<-function(sce,given_state){
    net<-sce@metadata$RL_net
    RL <- sce@metadata$RL_info
  Total_hover_list <- c("TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WW 0<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : WM 1<br>NPM1: MM 2<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WW 0<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: WM 1<br>NRAS: MM 2",
                        
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WW 0",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: WM 1",
                        "TET2 : MM 2<br>NPM1: MM 2<br>NRAS: MM 2")
  Total_label_list<-paste0("s",seq(0:(length(V(net))-1)))
  observed_states <-unique(RL$next_state[which(RL$observed_states==1)])
  clone_to_observed_clone <-igraph::shortest_paths(net,from=V(net)[given_state],to=V(net)[observed_states])$vpath
  V(net)[unlist(clone_to_observed_clone)]$color <-6
  g<- induced.subgraph(net,V(net)[unlist(clone_to_observed_clone)])
  visNet<-visIgraph(g)%>%
    #visNetwork(nodes=visNet$nodes,edges=visNet$edges,width="100%")%>%
    visNodes(color=ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey"),
             #label=paste0("s",seq(0:(length(V(net))-1)))
    )%>%
    visEdges(color="black")%>%
    visIgraphLayout(layout="layout_with_sugiyama")
  visNet$x$options$nodes$font=list(size=48)
  visNet$x$options$width="100%"
  visNet$x$nodes$label <-Total_label_list[unique(unlist(clone_to_observed_clone))]
  visNet$x$options$nodes$color<-ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  visNet$x$nodes$color<-ifelse(names(V(g))%in%droplevels(unique(RL$next_state[RL$observed_states==1])),"darkred","grey")
  
  visNet$x$nodes$title<-Total_hover_list[unique(unlist(clone_to_observed_clone))[order(unique(unlist(clone_to_observed_clone)),decreasing = FALSE)]]
  return(visNet)
}

