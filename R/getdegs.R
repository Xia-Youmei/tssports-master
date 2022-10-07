#' DESeq2 differential analysis
#' get differentially expressed genes (DEGs)
#' @param treat_group The 11th and 12th digits (last two digits) of the name of
#' the SRR treatment sample group were used to define the experimental
#' conditions and distinguish the control group from the experimental group.
#'
#' @return Form new file
#' "sports_counts_all.txt","sports_DEG_fdr005_2fc_all.txt",
#' "sports_cpm_fdr005_2fc_all.txt","sports_DEG_all.txt"
#' @export
#'
#' @examples getdegs(20:22)
getdegs <- function(treat_group) {
  library(tidyverse)
  suppressWarnings(suppressMessages(library(DESeq2)))
  counts=read.table("sports_combined_sample_fragments_counts_matrix_0.5.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
  for(i in 1:length(colnames(counts))){
  	colnames(counts)[i]=strsplit(colnames(counts)[i],"/")[[1]][1]
  }
  experiments_design=read.table("sports_combined_sample_fragments_annotation.txt",header = T, check.names = F)
  experiments_design <- experiments_design[!duplicated(experiments_design$Fragment),]
  rownames(experiments_design) <- experiments_design$Fragment
  experiments_design<-experiments_design[,-1]
  comgene <- intersect(rownames(counts),rownames(experiments_design))
  counts <- counts[comgene,]
  class(counts)
  class(comgene)
  experiments_design <- experiments_design[comgene,]
  a <- rownames(counts)
  b <- rownames(experiments_design)
  identical(a,b)
  counts$Gene <- as.character(experiments_design$Annotation)
  counts <- counts[!duplicated(counts$Gene),]
  rownames(counts) <- counts$Gene
  ncol(experiments_design)
  nrow
  counts <- counts[,-ncol(counts)]
  write.table(counts, file = "sports_counts_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  conditions=data.frame(sample=colnames(counts),
                        group=factor(ifelse(substr(colnames(counts),10,11) == c(treat_group),"T","N"),levels = c("N","T"))) %>%
    column_to_rownames("sample")
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = conditions,
    design = ~ group)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds)
  save(res,file = "LIHC_DEG.rda")
  res_deseq2 <- as.data.frame(res)%>%
    arrange(padj) %>%
    dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)
  res_deseq2=res_deseq2[!is.na(res_deseq2$padj),]
  res_deseq2$ID <- rownames(res_deseq2)
  res_deseq2 <- res_deseq2[,c(length(res_deseq2),1:(length(res_deseq2)-1))]
  write.table(res_deseq2,file = "sports_DEG_fdr005_2fc_all.txt",sep = "\t",row.names = FALSE,quote = F)
  dat <- log2(edgeR::cpm(counts)+1)
  dat <- as.data.frame(dat)
  dat$ID <- rownames(dat)
  dat <- dat[,c(length(res_deseq2),1:(length(res_deseq2)-1))]
  write.table(dat,file = "sports_cpm_fdr005_2fc_all.txt",sep = "\t",row.names = FALSE,quote = F)
  res_deseq2 <- as.data.frame(res)
  res_deseq2$ID <- rownames(res_deseq2)
  res_deseq2 <- res_deseq2[,c(length(res_deseq2),1:(length(res_deseq2)-1))]
  write.table(res_deseq2,file = "sports_DEG_all.txt",sep = "\t",row.names = FALSE,quote = F)
}
