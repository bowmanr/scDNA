#' Run Trajectory Analysis after extraction from SingleCellExperiment object
#'
#' @param sce is the SingleCellExperiment Object
#'
#' @return RL_analysis this returns a paths with metadata about the mutation order
#' @export
#' @examples
#' 
trajectory_analysis<-function(sce,use_ADO=TRUE){
  
  mutation_states<-length(unique(sce@metadata$Architecture$final_annot))
  # This builds the theoretical Markov decision process (MDP) of all possible mutation combinations
  print("Building MDP")
  
  adj_list<-BuildMDP(mutation_states,use_ADO=TRUE)
  # we attach the weights as rewards for our list to design the possible and likely MDP.
  print("Adding Weighted Edges")

  adj_list<-attach_weights(sce,adj_list)

  set.seed(281330800) #281-330-8004  
  print("Starting Optimization")
  RL_output <-mdp_Q_learning_with_linklist(adj_list,discount=0.9,N=10000)
  print("Finished Optimization")
  net <- igraph::graph_from_data_frame(data.frame(from=as.numeric(levels(RL_output$current_state))[RL_output$current_state],
                                          to=as.numeric(levels(RL_output$next_state))[RL_output$next_state],
                                          weight=abs(1/(RL_output$Q_values_normalized))),directed=TRUE)
  sce@metadata$RL_net <-net
  #order of single path from WildType to Dominant clone
  print("Getting Mutation Paths/Policies")
  print(class(RL_output))
  print(colnames(RL_output))
  WT_to_dom_clone <-igraph::shortest_paths(net,from="0",to=as.character(RL_output$next_state[which(RL_output$reward==max(RL_output$reward))[1]][1]),algorithm = "dijkstra")$vpath[[1]]
  WT_to_dom_clone_policy<-NULL
  for (iter in 1:(length(WT_to_dom_clone)-1)){
    check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(WT_to_dom_clone[[iter]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                              next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(WT_to_dom_clone[[iter+1]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
    action_idx<-order(as.data.frame(RL_output)%>%
                        dplyr::filter((current_state==WT_to_dom_clone[[iter]]$name) &(next_state==WT_to_dom_clone[[iter+1]]$name))%>%
                        dplyr::pull(Q_values_normalized),decreasing = TRUE)
    check_sub_var$reward<-(as.data.frame(RL_output)%>%
      dplyr::filter((current_state==WT_to_dom_clone[[iter]]$name) &(next_state==WT_to_dom_clone[[iter+1]]$name))%>%
      dplyr::pull(reward))[action_idx[1]]
    check_sub_var$mutation_taken<-(as.data.frame(RL_output)%>%
      dplyr::filter((current_state==WT_to_dom_clone[[iter]]$name) &(next_state==WT_to_dom_clone[[iter+1]]$name))%>%
      dplyr::pull(action_type))[action_idx[1]]
    
    WT_to_dom_clone_policy<-rbind(WT_to_dom_clone_policy,check_sub_var)
  }
  print("Done with 1")
  #order of single path from WildType to final theoretical clone
  WT_to_last_state <-igraph::shortest_paths(net,from="0",to=as.character(RL_output$current_state[dim(RL_output)[1]][1]),algorithm ="dijkstra")$vpath[[1]]
  WT_to_last_state_policy<-NULL
  for (iter in 1:(length(WT_to_last_state)-1)){
    check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(WT_to_last_state[[iter]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                              next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(WT_to_last_state[[iter+1]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
    action_idx<-order(as.data.frame(RL_output)%>%
      dplyr::filter((current_state==WT_to_last_state[[iter]]$name) &(next_state==WT_to_last_state[[iter+1]]$name))%>%
      dplyr::pull(Q_values_normalized),decreasing = TRUE)
    check_sub_var$mutation_taken<-(as.data.frame(RL_output)%>%
      dplyr::filter((current_state==WT_to_last_state[[iter]]$name) &(next_state==WT_to_last_state[[iter+1]]$name))%>%
      dplyr::pull(action_type))[action_idx[1]]
    check_sub_var$reward<-(as.data.frame(RL_output)%>%
      dplyr::filter((current_state==WT_to_last_state[[iter]]$name) &(next_state==WT_to_last_state[[iter+1]]$name))%>%
      dplyr::pull(reward))[action_idx[1]]
    WT_to_last_state_policy<-rbind(WT_to_last_state_policy,check_sub_var)
  }
print("Done with 2")
  # all wild types to dominant
  All_WT_to_dom<-igraph::all_simple_paths(net,from="0",to=as.character(RL_output$next_state[which(RL_output$reward==max(RL_output$reward))[1]][1]))
  All_WT_to_dom_policy =NULL
  for(big_list_iter in 1:length(All_WT_to_dom)){
    single_WT_to_dom_policy=NULL
    single_WT_to_dom <- All_WT_to_dom[[big_list_iter]]
    for (iter in 1:(length(single_WT_to_dom)-1)){

      check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(single_WT_to_dom[[iter]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                                next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(single_WT_to_dom[[iter+1]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
      action_idx<-order(as.data.frame(RL_output)%>%
                          dplyr::filter((current_state==single_WT_to_dom[[iter]]$name) &(next_state==single_WT_to_dom[[iter+1]]$name))%>%
                          dplyr::pull(Q_values_normalized),decreasing = TRUE)
      check_sub_var$mutation_taken<-(as.data.frame(RL_output)%>%
        dplyr::filter((current_state==single_WT_to_dom[[iter]]$name) &(next_state==single_WT_to_dom[[iter+1]]$name))%>%
        dplyr::pull(action_type))[action_idx[1]]
      check_sub_var$reward<-(as.data.frame(RL_output)%>%
        dplyr::filter((current_state==single_WT_to_dom[[iter]]$name) &(next_state==single_WT_to_dom[[iter+1]]$name))%>%
        dplyr::pull(reward))[action_idx[1]]
      single_WT_to_dom_policy<-rbind(single_WT_to_dom_policy,check_sub_var)
    }
    All_WT_to_dom_policy[[big_list_iter]]<-single_WT_to_dom_policy
  }
print("Done with 3")
  # set goal
  given_state="0"
  observed_states <-as.character(unique(RL_output$next_state[which(RL_output$observed_states==1)]))
  clone_to_observed_clone <-igraph::shortest_paths(net,from=V(net)[given_state],to=V(net)[observed_states])$vpath

  WT_to_all_observed_policy =NULL
  for(big_list_iter in 1:length(clone_to_observed_clone)){
    single_WT_to_obs_policy=NULL
    single_WT_to_obs <- clone_to_observed_clone[[big_list_iter]]
    if(length(single_WT_to_obs)>1){
      for (iter in 1:(length(single_WT_to_obs)-1)){
         check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(single_WT_to_obs[[iter]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                                   next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(single_WT_to_obs[[iter+1]]$name, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
         action_idx<-order(as.data.frame(RL_output)%>%
                             dplyr::filter((current_state==single_WT_to_obs[[iter]]$name) &(next_state==single_WT_to_obs[[iter+1]]$name))%>%
                             dplyr::pull(Q_values_normalized),decreasing = TRUE)
         check_sub_var$mutation_taken<-(as.data.frame(RL_output)%>%
          dplyr::filter((current_state==single_WT_to_obs[[iter]]$name) &(next_state==single_WT_to_obs[[iter+1]]$name))%>%
          dplyr::pull(action_type))[action_idx[1]]
         check_sub_var$reward<-(as.data.frame(RL_output)%>%
          dplyr::filter((current_state==single_WT_to_obs[[iter]]$name) &(next_state==single_WT_to_obs[[iter+1]]$name))%>%
          dplyr::pull(reward))[action_idx[1]]
        single_WT_to_obs_policy<-rbind(single_WT_to_obs_policy,check_sub_var)
      }
      WT_to_all_observed_policy[[big_list_iter]]<-single_WT_to_obs_policy
    } else {
      WT_to_all_observed_policy[[big_list_iter]]<-"No action taken: Starting clone is equal to ending clone"
    }
  }
  
  print("Saving Data")
  
  # need to reassign those from clone_code to mutation name order like Architecture id.
  sce@metadata$RL_info <-RL_output
  sce@metadata$Trajectories <- list("WT_to_dominant_clone"=WT_to_dom_clone_policy,
                                    "WT_to_last_possible_clone"=WT_to_last_state_policy,
                                    "All_WT_to_dominant_clone"= All_WT_to_dom_policy,
                                    "WT_to_Observed_paths"=WT_to_all_observed_policy)
  return(sce)
}
