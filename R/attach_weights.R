#' Assigns observed counts as weights to the Markov Decision Process
#'
#' @param sce The single cell experiment object that carries all data
#' @param adj_linklist this is an adjacency list for the theoretical states describing the MDP
#' @export
attach_weights <-function(sce,adj_linklist){
  
  # first we obtain weights from our final_summary dataframe (this was taken from the myeloid one)
  # TODO: Change the filtering, not sure if we always want "complete" or "other" or both?
  data_frame_clarity<-dplyr::inner_join(sce@metadata$Architecture,
                                       as.data.frame(sce@metadata$Clones)%>%
                                         dplyr::filter(Group=="Complete")%>%
                                                      dplyr::select(Clone,Count)%>%dplyr::distinct(),
                                      by="Clone")%>%
    dplyr::group_by(Clone)
  check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(adj_linklist$current_state, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                            next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(adj_linklist$next_state, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
  check2<-array2DF(t(rbind(mapply(function(x,y) unique(sce@metadata$Architecture$final_annot)[(which(bitwXor(utf8ToInt(x),utf8ToInt(y))>0)+1)/2], 
         check_sub_var$current_state,check_sub_var$next_state,SIMPLIFY = FALSE))))#%>%as_tibble()
  
  check_sub_var$mutation_taken<-"none"
  check_sub_var$mutation_taken[as.numeric(rownames(check2))]<-check2$Value
  
  weights_unique_list <-dplyr::distinct(data.frame(reward=data_frame_clarity$Count/sum(data_frame_clarity$Count)*100,
                                            names=as.numeric(gsub("_", "",data_frame_clarity$Clone))))
  
  colnames(weights_unique_list)<-c("reward","next_state")
  weights_unique_list$next_state<-as.factor(weights_unique_list$next_state)
  
  # here, we do three things:
  #   1) attach the weights as a reward,
  #   2) use these as observed clonograph states for searching later
  #   3) assign staying in place as a reward of 0
  adj_linklist <- dplyr::left_join(adj_linklist,weights_unique_list, by="next_state")%>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~tidyr::replace_na(.x,0)))
  
  # for search later
  adj_linklist$observed_states <- ifelse(adj_linklist$reward>0,1,0)
  adj_linklist$mutation_taken<-check_sub_var$mutation_taken
  # for ensuring we move forward through the MDP
  adj_linklist$reward[which(adj_linklist$current_state==adj_linklist$next_state)]=0
  adj_linklist$mutation_taken<-check_sub_var$mutation_taken
  
  adj_linklist<-adj_linklist%>%dplyr::select("current_state","mutation_taken","next_state","reward","Q_values","observed_states","legal_action","total_actions","Probability_transition")
  # Now we can run RL in the next function
  return(adj_linklist)
}

