#' A way to navigate the MDP for any given starting root node to leaf node.
#'
#'@param sce a singleCellExperiment object
#'@param start_name a clone_code of starting point. make sure it is less than goal_name
#'@param goal_name a clone_code for the end point to get to. make sure it is larger than goal_name
#'@export
#'@return It augments the sce object in the Trajectories list in the meta data by adding another path
get_own_path<-function(sce,start_name,goal_name){
  RL_output<-sce@metadata$RL_info
  # TODO: check if goal state is bigger than start state and swap if needed.
  start_name<-as.character(as.numeric(gsub("_", "",start_name)))
  goal_name<-as.character(as.numeric(gsub("_", "",goal_name)))
  clone_to_goal_clone <-igraph::shortest_paths(sce@metadata$RL_net,from=start_name,to=goal_name,algorithm = "dijkstra")$vpath[[1]]
  WT_to_state_policy<-NULL
  for (iter in 1:(length(clone_to_goal_clone)-1)){
    check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(levels(RL_output$next_state)[clone_to_goal_clone[iter]], length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                              next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(levels(RL_output$next_state)[clone_to_goal_clone[iter+1]], length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
    check_sub_var$reward<-as.data.frame(RL_output)%>%
      dplyr::filter((current_state==(clone_to_goal_clone)[iter][[1]]$name) &(next_state==(clone_to_goal_clone)[iter+1][[1]]$name))%>%
      dplyr::pull(reward)
    check_sub_var$mutation_taken<-as.data.frame(RL_output)%>%
      dplyr::filter((current_state==(clone_to_goal_clone)[iter][[1]]$name) &(next_state==(clone_to_goal_clone)[iter+1][[1]]$name))%>%
      dplyr::pull(action_type)
    WT_to_state_policy<-rbind(WT_to_state_policy,check_sub_var)
  }

  sce@metadata$Trajectories <-append(sce@metadata$Trajectories,
                                     list(WT_to_state_policy))
  names(sce@metadata$Trajectories)[length(sce@metadata$Trajectories)]<-paste0(start_name,"_to_",goal_name,"_path")
  return(sce)
}
