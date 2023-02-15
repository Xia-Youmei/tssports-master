#' @title Prepare tsRNA mapping plot data
#'
#' @description Prepare tsRNA mapping plot data, the data can use to plot the coverage profile of the tRNA.
#'
#'
#' @param tsRNAname  specific tsRNA name, which you can find in "_tRNA_mature_mapping.txt".
#' @param pattern  A uniform suffix for the filename results returned by the SPORTS (Shi et al.,2018) software.
#'
#' @return combined specific tsRNA mapping results, form a new file "specific tsRNA"_mappingplot.txt.
#' @export
#'
#' @examples prepare_tsrna_mapping_plotdata(tsRNAname = "Mus_musculus_mt_tRNA-Ala-TGC-1")
prepare_tsrna_mapping_plotdata<- function(tsRNAname,
                                           pattern = "_tRNA_mature_mapping.txt"){
  library(tidyverse)
  library(stringr)
  library(stringi)
  patternAddStar <- paste('*', pattern, '$', sep = '')
  filelist = list.files(pattern = patternAddStar)
  collapse_miRNA = lapply (filelist,function(x)read.table(x,sep = "\t", check.names = F,stringsAsFactors = F,header = F))
  cha1<-c(filelist)
  cha1 = str_split_fixed(cha1, pattern, 2)
  cha1 <- cha1[,1]
  aa <- list()
  aa3 <- list()
  tsRNAname <- tsRNAname
  tsRNAname1 <- paste(tsRNAname,'$',sep = '')
  for(i in 1:length(collapse_miRNA)){
    names(collapse_miRNA) <- cha1
    aa<-subset(collapse_miRNA[[i]],grepl(tsRNAname1,V1))
    if(tsRNAname %in% collapse_miRNA[[i]]$V1) aa2 <- c(tsRNAname,aa$V2,names(collapse_miRNA[i])) else aa2 <- c(tsRNAname,NA,names(collapse_miRNA[i]))
    aa3 <- rbind(aa3,aa2)
  }
  filename2 <- paste(tsRNAname,'_mappingplot.txt',sep = '')
  write.table(aa3,filename2,sep = "\t",row.names = FALSE,col.names=FALSE, quote = F)
}
