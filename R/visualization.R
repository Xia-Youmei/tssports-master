#' visualization_plot
#'
#' @return form 5 new file pie_plot_tsRNA_aa.pdf,pie_plot_tsRNA_end.pdf,maplot.pdf,heatmap_plot.pdf,volcano_plot.pdf.
#' @export
#'
#' @examples visualization()
visualization <- function(){
################################pie_plot_tsRNA_aa####################
  library(stringr)
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  tsRNA <- read.table("tsRNA_diff.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
  data<-data.frame(value=c(nrow(subset(tsRNA,grepl("Glu",tsRNA$ID))), nrow(subset(tsRNA,grepl("Gly",tsRNA$ID))), nrow(subset(tsRNA,grepl("Val",tsRNA$ID))), nrow(subset(tsRNA,grepl("Ser",tsRNA$ID)) )),group=c("Glu-tsRNA", "Gly-tsRNA", "Val-tsRNA", "Ser-tsRNA"))
  data_1 <- data %>% mutate(ratio = paste0(round(data$value/sum(data$value)*100,2),"%"))
  blank_theme <- theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position = 'right',
          legend.text = element_text(colour = 'black',size = 12),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 20)
    )
  p <- ggplot(data=data_1, mapping=aes(x="Improved",y=value,fill=group))+
    labs(x=NULL,y=NULL,title = 'N vs T')+
    geom_bar(stat="identity",width=1,color="white",position='stack',size=0.8)+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c( "#F58C63","#66C3A6","#8BA0CC","#A680B9"))+
    blank_theme +
    geom_text(stat="identity",aes(y=value, label = ratio), size=4, position=position_stack(vjust = 0.5))  #position_stack影响百分比在饼图内位置，加了才会正
  p
  ggsave(p,filename = 'pie_plot_tsRNA_aa.pdf',width=9,height=6)
  ################################pie_plot_tsRNA_end####################
  library(stringr)
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  tsRNA <- read.table("tsRNA_diff.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
  data<-data.frame(value=c(nrow(subset(tsRNA,grepl("5_end",tsRNA$ID))), nrow(subset(tsRNA,grepl("3_end",tsRNA$ID))), nrow(subset(tsRNA,grepl("CCA_end",tsRNA$ID)))),group=c("5'-tsRNA", "3'-tsRNA", "CCA-tsRNA"))
  data_1 <- data %>% mutate(ratio = paste0(round(data$value/sum(data$value)*100,2),"%"))
  blank_theme <- theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position = 'right',
          legend.text = element_text(colour = 'black',size = 12),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 20)
    )
  p <- ggplot(data=data_1, mapping=aes(x="Improved",y=value,fill=group))+
    labs(x=NULL,y=NULL,title = 'N vs T')+
    geom_bar(stat="identity",width=1,color="white",position='stack',size=0.8)+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c( "#F58C63","#66C3A6","#8BA0CC"))+
    blank_theme +
    geom_text(stat="identity",aes(y=value, label = ratio), size=4, position=position_stack(vjust = 0.5))  #position_stack影响百分比在饼图内位置，加了才会正
  p
  ggsave(p,filename = 'pie_plot_tsRNA_end.pdf',width=9,height=6)
  ################################maplot####################
  library(ggpubr)
  library(ggplot2)
  DEG <- read.table("sports_DEG_all.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
  diff_express <- DEG[c(3,1,2,6)]
  a <- rownames(diff_express)
  names(diff_express)=c("name","baseMean","log2FoldChange","padj")  # 新加一列Label，为了能在MA图中标记指定基因名
  diff_express$name <- a
  p <- ggmaplot(diff_express, main = expression(""), #expression("Group 1" %->% "Group 2")
                fdr = 0.05, fc = 2, size = 0.6,
                palette = c("#B31B21", "#1465AC", "darkgray"),
                genenames = as.vector(diff_express$name),
                legend = "top", top = 20,
                font.label = c("bold", 11), label.rectangle = TRUE,
                font.legend = "bold", select.top.method = "fc",
                font.main = "bold",
                ggtheme = ggplot2::theme_minimal())
  p
  ggsave("maplot.pdf",p,width=12,height=10)
  ################################heatmap_plot####################
  library(tidyverse)
  library(pheatmap)
  library(viridis)
  exp <- read.table("sports_cpm_fdr005_2fc_all.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
  exp <- as.matrix(exp)
  DEG <- read.table("sports_DEG_fdr005_2fc_all.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
  logFC_cutoff <- 1
  type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
  type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
  DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
  table(DEG$change)
  cg = rownames(DEG)[DEG$change !="NOT"]
  exp_diff <- exp[cg,]
  group_list=factor(ifelse(substr(colnames(exp),10,11) == c(20:22),"T","N"),levels = c("N","T"))
  annotation_col=data.frame(group=group_list)
  rownames(annotation_col)=colnames(exp_diff)
  p <- pheatmap(exp_diff,
                annotation_col=annotation_col,
                scale = "row",
                show_rownames = F,
                show_colnames =F,
                color=colorRampPalette(viridis(3))(1000),
                border_color=NA,
                #color = colorRampPalette(c("navy", "white", "red"))(50),
                cluster_cols =F,
                fontsize = 10,
                fontsize_row=3,
                fontsize_col=3)
  dev.off()
  ggsave(p,filename = 'heatmap_plot.pdf',width=9,height=6)
  ################################volcano_plot####################
  library(ggpubr)
  library(ggrepel)
  library(dplyr)
  d=read.table("sports_DEG_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
  d <- d[,c(3,7,1)]
  names(d)[1]=c("log2FC")
  names(d)[3]=c("label")
  d=d[!is.na(d[,2]),]
  d$log10=-log10(d[,2])  # col 2 fdr
  d$class="none"
  d[d[,2] <= 0.05 & d[,1] >= 1,]$class="UP"
  d[d[,2] <= 0.05 & d[,1] <= -1,]$class="DOWN"
  up_num=nrow(d[d$class == "UP",])
  down_num=nrow(d[d$class == "DOWN",])
  d$class <- as.factor(d$class) # col 3 class
  xval=ceiling(max(abs(d[,1])))
  data<- bind_rows(
    d %>%
      filter(class == 'UP') %>%
      arrange(padj,desc(abs(log2FC))) %>%
      head(5),
    d %>%
      filter(class == 'DOWN') %>%
      arrange(padj,desc(abs(log2FC))) %>%
      head(5)
  )
  colors <- c("UP"="#FC4E07", "none"="#8D8F8E", "DOWN"="#00AFBB")
  p=ggplot(data=d, aes(x=log2FC,y=log10,color=class)) +
    geom_point(data = d[d$class=="UP",],aes(y=log10,color="UP"))+
    geom_point(data = d[d$class=="none",],aes(y=log10,color="none"))+
    geom_point(data = d[d$class=="DOWN",],aes(y=log10,color="DOWN"))+
    scale_color_manual(values = colors)+
    scale_x_continuous(limits = c(-xval,xval),breaks=-xval:xval)+
    ylab(paste0("-log10 ",names(d)[2]))+xlab("Log2FoldChange")+
    theme_pubr(base_size = 12,border=TRUE)+geom_hline(yintercept=1.3, linetype="dashed",                                                                                        color = "black", size=0.5)+geom_point(data = d[d$label!="",])+
    annotate(geom = 'text', label = paste0('UP_Number: ', up_num), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+
    annotate(geom = 'text', label = paste0('DOWN_Number: ', down_num), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)+labs(color = "class")
  p
  p <- p+geom_label_repel(data = data,
                          aes(x=log2FC,y=log10,color=class,label=label),
                          fontface="bold",
                          color="black",
                          box.padding=unit(1, "lines"),
                          point.padding=unit(0.5, "lines"),
                          segment.colour = "black",segment.size = 0.5,segment.alpha = 0.5,max.overlaps = Inf)

  ggsave(paste0("volcano","_plot.pdf"),p,width=8,height=9)
}
