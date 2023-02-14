#' collapse_mature_mirna_reads
#'
#' @param x A folder address
#' It can extract -3p$|-5p$|miR characteristics of non-coding RNA,
#' sum it's reads and it will add to the original input file
#' to form new file with the '_output_collapse_miRNA.txt' at the end.
#' @return form new file with the '_output_collapse_miRNA.txt' at the end.
#' @export
#'
#' @examples collapse_mature_mirna_reads()
collapse_mature_mirna_reads <- function() {

  library(tidyverse)
  library(stringr)
  library(stringi)
  filelist = list.files(pattern = "*_output.txt$")
  collapse_miRNA = lapply (filelist,function(x)read.table(x,sep = "\t",check.names = F,stringsAsFactors = F,header = T))  #lapply函数 对列表、数据框数据集进行循环,输入为列表;function()指令来命名和创建函数
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
  #-----------------------------batch processing------------------------
  for(i in 1:length(collapse_miRNA)){
    notmml <-  subset(collapse_miRNA[[i]],!grepl("-miR-|-mir-|-let-",Annotation))
    mat <- "((-3p)$|(-5p)$|miR)"
    mi <- list()
    y=lapply(stri_split_regex(stri_reverse(collapse_miRNA[[i]]$Annotation), pattern = '[;\\s]+', n = 150), stri_reverse)
    for(j in 1:length(y)){
      mi[j] <- list(str_subset(y[[j]],mat))
    }
    mi <- mi[!duplicated(mi)]
    mi <- do.call(rbind,mi)
    mi <- c(mi[,1])
    mi <- mi[!duplicated(mi)]
    mi <- as.data.frame(mi)
    a <- list()
    collapse_miRNA35p <- list()
    for (j in 1:length(mi[,])) {
      a[j] <- list(stringr::str_which(collapse_miRNA[[i]]$Annotation,mi[j,]))
      x <- a[[j]]
      collapse_miRNA35p[j] <- list(sum(collapse_miRNA[[i]]$Reads[unlist(x)]))
    }
    collapse_miRNA35p <- do.call(rbind.data.frame, collapse_miRNA35p)
    collapse_miRNA35p <- data.frame(mi,collapse_miRNA35p)
    collapse_miRNA35p <- cbind( collapse_miRNA35p[1],collapse_miRNA35p[1],"NA", collapse_miRNA35p[2],"Yes", collapse_miRNA35p[1])
    cnames=c("ID","Sequence", "Length", "Reads", "Match_Genome", "Annotation")
    colnames(collapse_miRNA35p)=cnames
    collapse_miRNA[[i]] <- rbind(notmml,collapse_miRNA35p)
    filename2 <- paste('SRR',col11,'_output_collapse_miRNA.txt',sep = '')
    re=collapse_miRNA[[i]][,]
    write.table(re,filename2[i],sep = "\t",row.names = FALSE,quote = F)
  }
}
