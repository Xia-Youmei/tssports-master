#' PCA(principal component analysis)
#'
#' @param treat_group The 11th and 12th digits (last two digits) of the name of
#' the SRR treatment sample group were used to define the experimental
#' conditions and distinguish the control group from the experimental group.
#'
#' @return Form new file pca.pdf
#' @export
#'
#' @examples pca(20:22)
pca <- function(treat_group){
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  dm <- read.table("sports_counts_all.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
  sel = order(apply(dm, 1, var), decreasing=TRUE)[1:1000]
  dm2000=dm[sel,]
  data.pca <- PCA(as.data.frame(t(dm2000)), graph = FALSE)
  conditions=data.frame(sample=colnames(dm),
                        group=factor(ifelse(substr(colnames(dm),10,11) == c(treat_group),"T","N"),levels = c("N","T"))) %>%
    column_to_rownames("sample")
  groupinfo <- conditions
  groupinfo=groupinfo[names(dm),]
  groupinfo <- as.data.frame(groupinfo)
  p=fviz_pca_ind(data.pca,
                 geom.ind = "point",
                 col.ind = conditions$group,
                 addEllipses = TRUE,
                 ellipse.type = "confidence",
                 legend.title = "Groups"
  )+theme_bw(base_size = 12)+scale_shape_manual(values=seq(0,8))
  ggsave("pca.pdf",p,width=8,height=8)
}
