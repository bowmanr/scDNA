#' Title
#'
#' @param sce 
#' @param slot 
#'
#' @return
#' @export
#'
#' @examples
fcs_export<-function(sce,
                     slot="DSB_norm",
                     save_path="~/sample_test.fcs"){
  print(paste("Exporting protein data from", slot))
  
  if(slot=="DSB_norm"){
    data<-altExp(sce)@assays@data$DSB_norm
  } else if(slot=="CLR_norm"){
    data<-altExp(sce)@assays@data$CLR_norm
  } else if(slot=="raw"){
    data<-altExp(sce)@assays@data$Protein
  }

  dta<-t(data)%>%
      data.frame()%>%
      tibble::rownames_to_column(var="Cell")%>%
      inner_join(sce@metadata$NGT,by="Cell")%>%
      dplyr::mutate(Complete_genotype=ifelse(Group=="Complete",1,0))%>%
      dplyr::select(!c(Cell,Clone,Group))%>%
      as.matrix
  
  dta[,1:45] <-exp(1)^dta[,1:nrow(data)]
  
  mode(dta)<-"numeric"
  # you need to prepare some metadata
  meta <- data.frame(name=dimnames(dta)[[2]],
                     desc=paste(dimnames(dta)[[2]])
  )
  
  meta$range <- apply(apply(dta,2,range),2,diff)
  meta$minRange <- apply(dta,2,min)
  meta$maxRange <- apply(dta,2,max)
  
  # a flowFrame is the internal representation of a FCS file
  ff <- new("flowFrame",
            exprs=dta,
            parameters=AnnotatedDataFrame(meta))
  # now you can save it back to the filesystem
  write.FCS(ff,save_path)
}