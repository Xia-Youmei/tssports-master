#' combine_read_counts_from_sports
#'
#' @param x A folder address.
#'
#' @return form new file
#' sports_combined_sample_fragments_counts_matrix_all.txt
#' sports_combined_sample_fragments_counts_matrix_0.5.txt
#' sports_combined_sample_fragments_annotation.txt
#' @export
#'
#' @examples combine_read_counts()
combine_read_counts <- function() {

  library(tidyverse)
  library(stringr)
  library(stringi)
  library(dplyr)

  filelist = list.files(pattern = "*_miRNA.txt$")
  collapse_miRNA = lapply (filelist,function(x)read.table(x,sep = "\t",check.names = F,stringsAsFactors = F,header = T))
  filelist = list.files(pattern = "*.txt$")
  cha1<-c(filelist)
  col1<-str_extract_all(cha1,"\\d")
  i<-1
  while(i<=length(col1)){
    if(length(col1[[i]])==0) col1<-col1[-i] else i<-i+1
  }
  col11<-numeric(length(col1))
  for(i in 1:length(col1)){
    l1<-length(col1[[i]])
    l11<-c()
    for(j in 1:l1)
      l11<-paste(l11,col1[[i]][j],sep="")
    col11[i]<-as.numeric(l11)
  }
  col11<-col11[!duplicated(col11)]
  col11
  count_results <- list()
  reone <- list()
  for(i in 1:length(collapse_miRNA)){
    reone[[i]] <- collapse_miRNA[[i]][,c(2,4)]
    filename <- paste('SRR',col11,sep = '')
    cnames=c("Sequence", filename[i])  #列名
    colnames(reone[[i]])=cnames  #换列名
  }
  count_results <- reone %>% reduce(full_join, by = "Sequence")
  for (i in 1:ncol(count_results)){
    count_results[,i][is.na(count_results[,i])] <- 0
  }
  write.table(count_results,file = "sports_combined_sample_fragments_counts_matrix_all.txt",sep = "\t",row.names = FALSE,quote = F)
  count_results <- read.table("sports_combined_sample_fragments_counts_matrix_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
  count_results <- count_results[rowSums(count_results == 0) <= (ncol(count_results)-1)*0.5, ]
  write.table(count_results,file = "sports_combined_sample_fragments_counts_matrix_0.5.txt",sep = "\t",row.names = FALSE,quote = F)
  count_results_ann <- list()
  for(i in 1:length(collapse_miRNA)){
    count_results_ann[[i]] <- collapse_miRNA[[i]][,c(2,3,5,6)]
    cnames=c("Fragment","Length","Match_Genome","Annotation")
    colnames(count_results_ann[[i]])=cnames
  }
  count_results <- count_results_ann %>% reduce(full_join, by = c("Fragment","Length","Match_Genome","Annotation")) #列表内数据框一样的都合并，不一样的保留
  write.table(count_results,file = "sports_combined_sample_fragments_annotation.txt",sep = "\t",row.names = FALSE,quote = F)
}
