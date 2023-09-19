#' Title
#'
#' @param file path to h5 file
#' @param VAF_cutoff Fraction of cells that are mutated for initial filtering of variants (default 0.005, meaning 0.05%)
#' @param GT_cutoff Fraction of cells that are successfully genotyped for initial filtering (default 0.2, meaning 20%)
#' @param DP_cutoff minimum number of reads necessary for a reliable genotype call in a single cell (default: 10)
#' @param GQ_cutoff minimum genotype quality necessary for a reliable genotype call in a single cell (default: 30)
#' @param AF_cutoff Deviation from 0, 50, or 100% for a reliable call of WT, Het or Hom respectively (default: 25)
#'
#' @return list
#' @export
#'
#' @examples
loom_to_sce<-function(file,
                    GT_cutoff=35,
                    VAF_cutoff=5,
                    DP_cutoff=10,
                    GQ_cutoff=30,
                    AF_cutoff=25,
                    variant_set=NULL){
  
file  
print(paste("Input file:",file))
print("Checking for Barcode Duplicates")
all_barcodes<- rhdf5::h5read(file=file,name="/col_attrs/barcode")
dedup_barcodes<-all_barcodes[!(duplicated(all_barcodes) | duplicated(all_barcodes, fromLast = TRUE))]
viable_barcodes<-which(rhdf5::h5read(file=file,name="/col_attrs/barcode")%in%dedup_barcodes)
print(paste(length(all_barcodes)-length(viable_barcodes), "duplicated barcodes detected and removed"))

if(is.null(variant_set)){
  print("Identifying Variants Past Threshold")
  print(paste("Genotyping Threshold",GT_cutoff))
  print(paste("VAF Threshold",VAF_cutoff, "For Any Sample (if merged)"))
  
    print(paste("Loading Genotype Data for:",sample))

    NGT_data<-rhdf5::h5read(file=file,name="/matrix")%>%
      {as(.,"dgCMatrix")}%>%
      `rownames<-`(., rhdf5::h5read(file=file,name="/col_attrs/barcode"))  %>% 
      `colnames<-`(., rhdf5::h5read(file=file,name="/row_attrs/id"))
    
    VAF_cut_variants<-data.frame("id" = rhdf5::h5read(file=file,name="/row_attrs/id"), 
                      "WT"=sparseMatrixStats::colCounts(NGT_data,value = 0),
                      "Het"=sparseMatrixStats::colCounts(NGT_data,value = 1),
                      "Hom"=sparseMatrixStats::colCounts(NGT_data,value = 2),
                      "Missing"=sparseMatrixStats::colCounts(NGT_data,value = 3))%>%
             dplyr::mutate(VAF=((Het+Hom*2)/((WT+Het+Hom)*2))*100)%>%
             dplyr::mutate(genotyping_rate=((WT+Het+Hom)/(WT+Het+Hom+Missing) )*100)%>%
             dplyr::filter(genotyping_rate>=GT_cutoff)%>%
             dplyr::filter(VAF>=VAF_cutoff)%>%
             dplyr::arrange(desc(genotyping_rate))
    
    
  if(return_variants_only==TRUE){
    return(VAF_cut_variants)
    break
  } 
  VAF_cut_variants_complete<-do.call(rbind,VAF_cut_variants)
  VAF_cut_variants<-VAF_cut_variants_complete%>%dplyr::pull(id)
  print(paste("Loading n=",nrow(VAF_cut_variants_complete),"variants")) 
  print(paste("VAF Range:",paste(round(range(VAF_cut_variants_complete$VAF),digits = 3) ,sep=" ",collapse = " ")))
}



if(!is.null(variant_set)){
  VAF_cut_variants<- variant_set
  print(paste("Loading n=",length(VAF_cut_variants),"variants")) 
  print(paste("VAF Range: Not currently computed for predefined variants"))
}

VAF_cut_index<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id")%in%VAF_cut_variants)
print(VAF_cut_index)
print("Loading Allele Frequency Data")           
AF_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/AF",index=list(VAF_cut_index,viable_barcodes))%>%
  `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
  data.frame()%>%
  dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
  tidyr::pivot_longer(cols=!id,values_to = "AF",names_to = "barcode")

print("Loading Depth Data")
print(paste("Depth Cutoff >",DP_cutoff))
DP_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/DP",index=list(VAF_cut_index,viable_barcodes))%>%
  `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
  data.frame()%>%
  dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
  tidyr::pivot_longer(cols=!id,values_to = "DP",names_to = "barcode")


print("Loading Genotype Quality Data")
print(paste("Genotype quality cutoff >",GQ_cutoff))
GQ_data<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/GQ",index=list(VAF_cut_index,viable_barcodes))%>%
  `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
  data.frame()%>%
  dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
  tidyr::pivot_longer(cols=!id,values_to = "GQ",names_to = "barcode")

print("Loading Subsetted Genotype Information")
NGT_data_subset<-rhdf5::h5read(file=file,name="/assays/dna_variants/layers/NGT",index=list(VAF_cut_index,viable_barcodes))%>%
  `colnames<-`(., rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)))  %>% 
  data.frame()%>%
  dplyr::mutate(id=rhdf5::h5read(file=file,name="/assays/dna_variants/ca/id",index=list(VAF_cut_index))) %>%
  tidyr::pivot_longer(cols=!id,values_to = "NGT",names_to = "barcode")

print("Final Filtering")
sample_names<-data.frame("barcode"=rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode",index=list(viable_barcodes)),    
                         "Sample"=rhdf5::h5read(file=file,name="/assays/dna_variants/ra/sample_name",index=list(viable_barcodes)))%>%
  dplyr::mutate(barcode=gsub("-","\\.",barcode))


out<-list(AF_data,
          GQ_data,
          DP_data,
          NGT_data_subset)%>%
  purrr::reduce(dplyr::full_join, by=c("id","barcode"))%>%
  dplyr::full_join(sample_names,by="barcode")%>%
  dplyr::mutate(NGT_mask=ifelse(NGT==3,FALSE,TRUE))%>%
  dplyr::mutate(AF_mask=case_when( 
    (NGT==0 & AF<(50-AF_cutoff))~TRUE,
    (NGT==1 & (AF>(50-AF_cutoff)|AF<(50+AF_cutoff)))~TRUE,
    (NGT==2 & AF>(50+AF_cutoff))~TRUE,
    TRUE~FALSE
  ))%>%
  dplyr::mutate(DP_mask=ifelse(DP>DP_cutoff,TRUE,FALSE))%>%
  dplyr::mutate(GQ_mask=ifelse(GQ>GQ_cutoff,TRUE,FALSE))


final_barcodes<-out$barcode%>%unique%>%{gsub("\\.","-",.)}
final_barcodes_index<-which(rhdf5::h5read(file=file,name="/assays/dna_variants/ra/barcode")%in%final_barcodes)

print("Formating Single Cell Experiment")


sce<-SingleCellExperiment::SingleCellExperiment(list(
  NGT=tidyr::pivot_wider(out,
                         id_cols=id,
                         names_from=barcode,
                         values_fill=NA,
                         values_from = NGT)%>%
    dplyr::select(-id),
  NGT_mask=tidyr::pivot_wider(out,
                              id_cols=id,
                              names_from=barcode,
                              values_fill=NA,
                              values_from = NGT_mask)%>%
    dplyr::select(-id),                              
  AF=tidyr::pivot_wider(out,
                        id_cols=id,
                        names_from=barcode,
                        values_fill=NA,
                        values_from = AF)%>%
    dplyr::select(-id),
  AF_mask=tidyr::pivot_wider(out,
                             id_cols=id,
                             names_from=barcode,
                             values_fill=NA,
                             values_from = AF_mask)%>%
    dplyr::select(-id),
  DP=tidyr::pivot_wider(out,
                        id_cols=id,
                        names_from=barcode,
                        values_fill=NA,
                        values_from = DP)%>%
    dplyr::select(-id),
  DP_mask=tidyr::pivot_wider(out,
                             id_cols=id,
                             names_from=barcode,
                             values_fill=NA,
                             values_from = DP_mask)%>%
    dplyr::select(-id),
  GQ=tidyr::pivot_wider(out,
                        id_cols=id,
                        names_from=barcode,
                        values_fill=NA,
                        values_from = GQ)%>%
    dplyr::select(-id),
  GQ_mask=tidyr::pivot_wider(out,
                             id_cols=id,
                             names_from=barcode,
                             values_fill=NA,
                             values_from = GQ_mask)%>%
    dplyr::select(-id)
))

SummarizedExperiment::colData(sce)<-S4Vectors::DataFrame(sample=out%>%
                                                           dplyr::distinct(barcode,Sample)%>%
                                                           pull(Sample))
rownames(sce)<-tidyr::pivot_wider(out,
                                  id_cols=id,
                                  names_from=barcode,
                                  values_fill=NA,
                                  values_from = NGT)%>%
  dplyr::pull(id)
colnames(sce)<-tidyr::pivot_wider(out,
                                  id_cols=id,
                                  names_from=barcode,
                                  values_fill=NA,
                                  values_from = NGT)%>%
  dplyr::select(-id)%>%
  colnames()
print("Adding on protein Data")
if(protein==TRUE){
  protein_mat <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/layers/read_counts",index=list(NULL,viable_barcodes))
  rownames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ca/id")
  colnames(protein_mat) <- rhdf5::h5read(file = file, name = "/assays/protein_read_counts/ra/barcode",index=list(viable_barcodes))
  SingleCellExperiment::altExp(sce, "Protein") <- SingleCellExperiment::SingleCellExperiment(list(Protein=protein_mat))
}
return(sce)
}

