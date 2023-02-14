#' tsRNA_ann_classify
#'
#' @return form new file
#' "miRNA_diff.txt","tsRNA_diff.txt",
#' "rsRNA_diff.txt","ysRNA_diff.txt"
#' @export
#'
#' @examples tsRNA_ann_classify()
tsRNA_ann_classify <- function(){
library(tidyverse)
library(stringr)
library(stringi)

collapse_miRNA <- read.table("sports_DEG_fdr005_2fc_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
miRlet <-  subset(collapse_miRNA,grepl("-miR-|-let-",collapse_miRNA$ID))
tsRNA <-  subset(collapse_miRNA,grepl("tRNA",collapse_miRNA$ID))
rsRNA <-  subset(collapse_miRNA,grepl("rRNA",collapse_miRNA$ID))
ysRNA <-  subset(collapse_miRNA,grepl("YRNA",collapse_miRNA$ID))
write.table(miRlet,"miRNA_diff.txt",sep = "\t",row.names = FALSE,quote = F)
write.table(tsRNA,"tsRNA_diff.txt",sep = "\t",row.names = FALSE,quote = F)
write.table(rsRNA,"rsRNA_diff.txt",sep = "\t",row.names = FALSE,quote = F)
write.table(ysRNA,"ysRNA_diff.txt",sep = "\t",row.names = FALSE,quote = F)
}
