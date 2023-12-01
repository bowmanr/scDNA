#' cell_confidence_labeling 
#'  This function provides additional labeling to he quality of the cells based 
#'  on dna and protein data. 
#' @param sce object 
#'
#' @return sce object with updated labeling in the colData to be used later
#' @export
#'
#' @examples
cell_confidence_labeling<-function(sce){


all_protein_droplets <- rhdf5::h5read(file = sce@metadata$file, 
                                      name = "/all_barcodes/protein_read_counts/layers/read_counts")
all_dna_droplets <- rhdf5::h5read(file = sce@metadata$file, 
                                  name = "/all_barcodes/dna_read_counts/layers/read_counts")
colnames(all_dna_droplets) <- rhdf5::h5read(file = sce@metadata$file, 
                                            name = "/all_barcodes/dna_read_counts/ra/barcode")
colnames(all_protein_droplets) <- rhdf5::h5read(file = sce@metadata$file, 
                                                name = "/all_barcodes/protein_read_counts/ra/barcode")
rownames(all_protein_droplets) <-rhdf5::h5read(file=sce@metadata$file,
                                               name="/all_barcodes/protein_read_counts/ca/id")

all_dna <- data.frame(Cell=colnames(all_dna_droplets),
                      Amps=t(all_dna_droplets))

all_protein<-data.frame(Cell=colnames(all_protein_droplets),
                        t(all_protein_droplets))

df2 <- dplyr::inner_join(all_dna, all_protein)
df2<-df2%>%dplyr::mutate(Cell = gsub("-1","", Cell))
df2<-df2 %>% dplyr::mutate(AMP_tot=rowSums(df2[,2:dim(all_dna)[2]]),
                    AMP_hits = rowSums(df2[,2:dim(all_dna)[2]]>0),
                    Pro_tot=rowSums(df2[,(dim(all_dna)[2]+1):(dim(all_dna)[2]+dim(all_protein)[2]-1)]),
                    Pro_hit=rowSums(df2[,(dim(all_dna)[2]+1):(dim(all_dna)[2]+dim(all_protein)[2]-1)]>0))
                    
raw_data_to_check = cbind((df2$AMP_tot),df2$AMP_hits,(df2$Pro_tot),df2$Pro_hit)
# use kernel density as our whack guassian mixture model
tab<-table(df2$AMP_hits)
cut_x=order(head(tab,-round(0.1*(length(tab)))),decreasing = FALSE)[1]

#print(paste("this is total cell count: ",dim(raw_data_to_check)[1]))
#print(paste("this is minimum amplicons needed: ",cut_x))
rownames(df2)<-df2$Cell
# Make cut to minimum point. DO I include 0s in left_cut? hard to say right now.
right_cut=raw_data_to_check[raw_data_to_check[,2]>=cut_x,]
right_cut_cells = df2[raw_data_to_check[,2]>=cut_x,]

test_point = rep(1,dim(all_dna)[2]-1) 
test_prob = test_point/sum(test_point)
right_cut_amplicon_counts<-(right_cut_cells[,2:(dim(all_dna)[2])])

right_cut_amplicon_dist <- right_cut_amplicon_counts/rowSums(right_cut_amplicon_counts)
Entropy_amplicons <- (-1)*rowSums(right_cut_amplicon_dist*log2(right_cut_amplicon_dist),na.rm=TRUE)
Entropy_max = -sum(test_prob*log2(test_prob),na.rm=TRUE)
Efficiency = Entropy_amplicons/Entropy_max
new_hueristic2 = rowSums(right_cut_amplicon_counts)*Efficiency
right_cut_cells$score <-new_hueristic2
dim(right_cut_cells)

reduced_df <-right_cut_cells%>%
  dplyr::filter(Cell%in%sce@colData@rownames)

# protein entropy!
protein_mat <- reduced_df[,(dim(all_dna)[2]+1):(dim(all_dna)[2]+dim(all_protein)[2]-1)]

protein_prob <- log10(protein_mat+1)/rowSums(log10(protein_mat+1))

protein_entropy <- (-1)*rowSums(protein_prob*log2(protein_prob),na.rm=TRUE)
protein_prob
prob_sorted = t(apply(protein_prob,1,sort))
prob_sorted
mySum_sorted = t(apply(prob_sorted, 1, cumsum))

uniform_list <- cumsum(rep(1/(dim(mySum_sorted)[2]),(dim(mySum_sorted)[2])))
uniform_diff<-rbind(rep(1/(dim(mySum_sorted)[2]),(dim(mySum_sorted)[2]-1)),diff(uniform_list))
uniform_unit<-apply(uniform_diff,2,function(x){
  y = c(1,0)
  dot.prod <- x%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  heading<-acos(dot.prod/(norm.x*norm.y))
  return(heading)})

xte <- matrix(NA,nrow=dim(mySum_sorted)[1],ncol=dim(mySum_sorted)[2])
cell_heading <- matrix(NA,nrow=dim(mySum_sorted)[1],ncol=(dim(mySum_sorted)[2]-1))
for (iter in 1:dim(mySum_sorted)[1]){
  cell_positions <-rbind(uniform_list,(mySum_sorted[iter,]))
  for(iter_j in 1:dim(cell_positions)[2]){
    doublet_point <-rbind(uniform_list[iter_j],uniform_list[iter_j])
    xte[iter,iter_j]<-sqrt(sum((cell_positions[,iter_j]-doublet_point)**2))
  }
  cell_trajectory <- rbind(rep(1/(dim(mySum_sorted)[2]),(dim(mySum_sorted)[2]-1)),diff(mySum_sorted[iter,]))
  cell_heading[iter,] <-apply(cell_trajectory,2,function(x){
    y = c(1,0)
    dot.prod <- x%*%y
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    heading<-acos(dot.prod/(norm.x*norm.y))
    return(heading)})
}
left_ecdf_func <- function(data) { 
  Length <- length(data) 
  sorted <- sort(data) 
  ecdf <- rep(0, Length) 
  for (i in 1:Length) { 
    ecdf[i] <- sum(sorted <= data[i]) / Length 
  } 
  return(ecdf) 
} 
right_ecdf_func <- function(data) { 
  Length <- length(data) 
  sorted <- sort(data) 
  ecdf <- rep(0, Length) 
  for (i in 1:Length) { 
    ecdf[i] <- sum(sorted >= data[i]) / Length 
  } 
  return(ecdf) 
} 

measures<-cbind(xte,cell_heading,reduced_df$score,reduced_df$Pro_tot)

outlier_score<-matrix(NA,nrow=dim(measures)[1],ncol = dim(measures)[2])
skew_meas <- apply(measures,2,skewness)
skew_meas[is.nan(skew_meas)]<-0
for(pro_iter in 1:dim(measures)[2]){
  right_ecdf<-right_ecdf_func(measures[,pro_iter])
  left_ecdf<-left_ecdf_func(measures[,pro_iter])
  for(cell_iter in 1:dim(measures)[1]){
    
    if(skew_meas[pro_iter]>0){
      outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
    }
    else{
      outlier_score[cell_iter,pro_iter]<-left_ecdf[cell_iter]
    }
  }
}
outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))

outlier_cutoff_val<-quantile(outlier_score_neglog)[4]+1.5*IQR(outlier_score_neglog)

reduced_df$outlier_score<-outlier_score_neglog

reduced_df$ent_score<-protein_entropy
colData(sce)$dna_score <-reduced_df$score
colData(sce)$protein_entropy<-reduced_df$ent_score
colData(sce)$outlier_score <-reduced_df$outlier_score
colData(sce)$cell_label_confidence <- ifelse(reduced_df$outlier_score<=outlier_cutoff_val,"High_confidence","Low_confidence")
colData(sce)$cell_label_confidence<-ifelse(reduced_df$score>=mean(new_hueristic2),colData(sce)$cell_label_confidence,"Poor_DNA")

return(sce)



}