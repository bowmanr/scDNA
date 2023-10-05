#' The main function builds an adjacency list of the theoretical order of mutations. The same MDP mask can be reused different mutations
#' 
#' This contains 1 main function and 3 helper functions.
#' The 3 helper functions aim to:
#' 1) turn our states to Ternary variables (can be 0(WT), 1(Het), or 2(Hom))
#' 2) a function that finds the difference in states that are < 1 mutation away
#' 3) a function that does bit logic subtraction for the ternary variables
#' @export
#' @param num_mutations The number of variants we are using (automatically obtained from the Architecture)

BuildMDP <-function(num_mutations,use_ADO=FALSE){

  num_mutations =3
  getidx <-function(n) {ifelse( (n==0|n==1),TRUE,FALSE)}
  dsum <- function(n) {ifelse(n < 10, n, n %% 10 + dsum(floor(n / 10)))}
  getTernaryFromState<-function(x,base=3){
    ifelse(x<base,x,getTernaryFromState(x%/% base)*10+x%%base)
  }
  backwardCheck<-function(x){(ifelse(dsum((-1)*x)==1,TRUE,FALSE))
  }

  
  if(num_mutations>=9){
    load(system.file(paste0('data/AdjacencyList_',num_mutations,'mutations.rDa'), package = 'scDNA'))
    return(adj_list)
  }
  num_type = 3
  num_states =num_type^num_mutations
  # this builds our masking function, we loop over upper triangle.
  # check if the summed digits are different by 1, if so it is an ok state we accept
  # EXAMPLE, check if StateA subtracted by stateC and stateB
  # stateC 121, stateB 110
  # stateA 010
  # subtract = 111 , 100
  # now sum digits 1+1+1 = 3, so reject, vs. 1+0+0 = 1 so accept!
  temp_vals <-(data.frame(x=(0:(num_states-1)))%>%apply(.,1,getTernaryFromState))
  
  j<-NULL
  i<-NULL
  state_tern = temp_vals
  
  f.compfun3<-function(state_tern,num_muts){
    f<-compiler::cmpfun(function(state_tern=temp_vals,num_muts=num_mutations){
      j<-NULL
      i<-NULL
      state_tern = temp_vals
      print(length(state_tern))
      r<-foreach::foreach(state_to_check=1:length(state_tern), .combine='cbind', .multicombine=TRUE,.export=c("dsum","getidx")) %dopar%{
        temp<-matrix(NA_real_,nrow=2)
        state_val <-(dsum(state_tern-state_tern[state_to_check]))
        vec<-which(getidx(state_val))
        ADO <-which(ifelse(dsum((-1)*state_val)==1,TRUE,FALSE))
        check1 <- unlist(lapply(ADO,
                                function(x) ifelse((intToUtf8(
                                  utf8ToInt(
                                    stringr::str_pad(state_tern[state_to_check], 
                                                     num_muts, 
                                                     pad = "0"))[(which((utf8ToInt(stringr::str_pad(state_tern[x], 
                                                                                                    num_muts, 
                                                                                                    pad = "0"))-utf8ToInt(stringr::str_pad(state_tern[state_to_check], 
                                                                                                                                           num_muts, 
                                                                                                                                           pad = "0"))) !=0))]))=="1",x,NA)))
        forwardADO <- unlist(lapply(vec,
                                    function(x) ifelse((intToUtf8(
                                      utf8ToInt(
                                        stringr::str_pad(state_tern[state_to_check], 
                                                         num_muts, 
                                                         pad = "0"))[(which((utf8ToInt(stringr::str_pad(state_tern[x], 
                                                                                                        num_muts, 
                                                                                                        pad = "0"))-utf8ToInt(stringr::str_pad(state_tern[state_to_check], 
                                                                                                                                               num_muts, 
                                                                                                                                               pad = "0"))) !=0))]))=="1",x,NA)))
        
        
        vec<-append(vec,check1[!is.na(check1)])
        vec<-append(vec,forwardADO[!is.na(forwardADO)])
        print(state_to_check)
        j<-t(vec)
        i<-t(rep(state_to_check,length(vec)))
        temp <-rbind(i,j)
        return(temp)
      }
      r
    })
    f(f(0))
  }
  
  f.compfun2<-function(state_tern){
    f<-compiler::cmpfun(function(state_tern=temp_vals){
      j<-NULL
      i<-NULL
      state_tern = temp_vals
      print(length(state_tern))
      r<-foreach::foreach(state_to_check=1:length(state_tern), .combine='cbind', .multicombine=TRUE,.export=c("dsum","getidx")) %dopar%{
        temp<-matrix(NA_real_,nrow=2)
        vec<-which(getidx(dsum(state_tern-state_tern[state_to_check])))
        print(state_to_check)
        j<-t(vec)
        i<-t(rep(state_to_check,length(vec)))
        temp <-rbind(i,j)
        return(temp)
      }
      r
    })
    f(f(0))
  }
 
  # some parallelization clusters, and records time.
  if(use_ADO){
    cl <- parallel::makeCluster(7)
    doParallel::registerDoParallel(cl)
    #start<-Sys.time()
    output_mat <-f.compfun3(temp_vals,num_mutations)
    #print( Sys.time() - start )
    parallel::stopCluster(cl)
  }else{
    cl <- parallel::makeCluster(7)
    doParallel::registerDoParallel(cl)
    #start<-Sys.time()
    output_mat <-f.compfun2(temp_vals)
    #print( Sys.time() - start )
    parallel::stopCluster(cl)
  }

  # Now we have a 2 row matrix we want to make a sparse representation
  reward_mat<-Matrix::sparseMatrix(output_mat[1,],output_mat[2,],,1)
  
  reward_mat
  rownames(reward_mat) <-temp_vals
  colnames(reward_mat) <-temp_vals
  # this next line is basically a pivot longer but ditches all 0s
  adj_list<-mefa4::Melt(reward_mat)
  colnames(adj_list) <-c("current_state","next_state","legal_action")
  # this unsets the duplicate rows so we can set forward ADO and mutation
  adj_list<-adj_list%>%
    uncount(legal_action)%>%
    mutate(legal_action=1)
  
  # labeling forward ADO, backward ADO, mutation, or none
  adj_list$action_type <-"mutation"
  adj_list$action_type <-ifelse(dsum(as.numeric(adj_list$next_state)-as.numeric(adj_list$current_state))<0,"ADO","mutation")
  adj_list$action_type <-ifelse(dsum(as.numeric(adj_list$next_state)-as.numeric(adj_list$current_state))==0,"none",adj_list$action_type)
  adj_list<-adj_list %>%
    dplyr::group_by(current_state,next_state,legal_action,action_type)%>%
    dplyr::mutate(rank = rank(action_type,ties.method="first"))%>%
    dplyr::mutate(action_type=ifelse(rank==2,"forward_ADO",action_type))#%>%
    #dplyr::ungroup()%>%
    #dplyr::select(-rank)
    
    
    
  # For now we will say equal transition probability
  adj_list%>%
    dplyr::group_by(current_state)%>%
    dplyr::mutate(total_actions=sum(legal_action))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(Probability_transition=legal_action/total_actions)->adj_list
  
  #assign Q-values to be 0 for each one right now, these get updated in the mdp_Q_learning_with_linklist.R function
  adj_list$Q_values<-rep(0,length(adj_list$current_state))

  # This is where we would save our different sized linked lists as .rDa files so we can
  #mkdir(paste0(getwd(),"/Results"))
  #save(adj_list,file=paste0("./Results/AdjacencyList_",num_mutations,"_mutations.rDa"))
  return(adj_list)
}
