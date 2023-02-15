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

#' @title Plot tsRNA mapping comparision coverage profile
#'
#' @description
#'  Plot the coverage profile of the tsRNA, Solid curves presenting average RPM of 2 groups(treat VS control) or 3 groups(high VS control; low VS control). RPM, reads per million; nt, nucleotide. Colored band stood for 95% CI.
#'
#' @param filename  specific tsRNA name, which filename maybe like "_mappingplot.txt".
#' @param pattern  A uniform suffix for the filename results returned by the function of prepare_tsRNA_mapping_plotdata.
#' @param Tsamples A collection of treat group file name prefixes, generally after removing "_mappingplot.txt".
#' @param Hsamples A collection of high group file name prefixes, generally after removing "_mappingplot.txt".
#' @param Lsamples A collection of low group file name prefixes, generally after removing "_mappingplot.txt".
#' @param Nsamples A collection of control group file name prefixes, generally after removing "_mappingplot.txt".
#' @param colors The color of the line in the line chart.
#' @param suffix A format suffix for a graph file.
#'
#' @return the coverage profile of the tsRNA and form a new file "specific tsRNA".pdf.
#' @export
#'
#' @examples
#' 2 group : treat VS control
#' tsRNA_mapping_comparision(filename = "Homo_sapiens_tRNA-Glu-TTC-1_mappingplot.txt",
#' Tsamples=c("XWH5_S7_L004", "XWH4_S6_L004", "XWH1_S0_L000", "XWH2_S6_L003"),
#' Nsamples=c("XWN1_S1_L004", "XWN2_S2_L004", "XWN3_S5_L003")
#'
#' 3 group : high VS control; low VS control
#' tsRNA_mapping_comparision(filename = "Homo_sapiens_tRNA-Glu-TTC-1_mappingplot.txt",
#' Hsamples=c("XWH5_S7_L004", "XWH4_S6_L004", "XWH1_S0_L000", "XWH2_S6_L003"),
#' Lsamples=c("XWH3_S5_L004",  "XWL4_S11_L004", "XWL1_S8_L004",  "XWL3_S0_L000",  "XWL2_S9_L004",  "XWL5_S1_L001" ),
#' Nsamples=c("XWN1_S1_L004", "XWN2_S2_L004", "XWN3_S5_L003"))
#'
#'
#'
tsRNA_mapping_comparision<- function(filename,
                                     Tsamples= FALSE,
                                     Hsamples= FALSE,
                                     Lsamples= FALSE,
                                     Nsamples,
                                     pattern = "_mappingplot.txt",
                                     colors = c("#00AFBB","#FC4E07","#E7B800"),
                                     suffix='.pdf'){
  library(wrMisc)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)

  tRNA.mapping <- read.table(filename ,sep="\t")
  tRNA.number <- nrow(tRNA.mapping)

  tsRNAname = str_split_fixed(filename, pattern, 2)
  tsRNAname <- tsRNAname[,1]

  re.data=data.frame(length=1:length(strsplit(as.character(tRNA.mapping[1,2]), ",")[[1]]))
  for (i in 1:tRNA.number) {
    rpm <- strsplit(as.character(tRNA.mapping[i,2]), ",")
    title.name <- as.character(tRNA.mapping[i,1])
    len <- length(rpm[[1]])
    length = c(1:len)
    length.data <- data.frame(RPM = as.numeric(rpm[[1]]))
    colnames(length.data)[1]=tRNA.mapping[i,3]
    re.data=cbind(re.data,length.data)
  }

  if( Tsamples != FALSE ){
    Tsamples=Tsamples
    Nsamples=Nsamples

    tsem=rowSEMs(re.data[,Tsamples])
    nsem=rowSEMs(re.data[,Nsamples])

    tmean=rowMeans(re.data[,Tsamples])
    nmean=rowMeans(re.data[,Nsamples])

    re.data=data.frame(length=length,tmean=tmean,nmean=nmean)

    dp=melt(re.data, id.vars = "length")

    dp$sem=c(tsem,lsem,nsem)
    dp$low=dp$value-dp$sem
    dp$high=dp$value+dp$sem
    dp$variable=factor(dp$variable,levels=c("nmean","tmean"))
  }  else {
    Hsamples=Hsamples
    Lsamples=Lsamples
    Nsamples=Nsamples

    hsem=rowSEMs(re.data[,Hsamples])
    lsem=rowSEMs(re.data[,Lsamples])
    nsem=rowSEMs(re.data[,Nsamples])

    hmean=rowMeans(re.data[,Hsamples])
    lmean=rowMeans(re.data[,Lsamples])
    nmean=rowMeans(re.data[,Nsamples])

    re.data=data.frame(length=length,hmean=hmean,lmean=lmean,nmean=nmean)

    dp=melt(re.data, id.vars = "length")

    dp$sem=c(hsem,lsem,nsem)
    dp$low=dp$value-dp$sem
    dp$high=dp$value+dp$sem
    dp$variable=factor(dp$variable,levels=c("nmean","hmean","lmean"))
  }

  dp<-dp[complete.cases(dp),]
  p=ggplot(dp, aes(x = length, y = value, color=variable))+
    geom_line(size=1)+theme_pubr()+
    geom_ribbon(aes(ymin = low, ymax = high,fill=variable), alpha = 0.1,colour = NA)+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    ylab("Coverage(RPM)")+xlab("Length of nucleotides")+labs(color="group")+labs(fill="group")
  ggsave(paste(tsRNAname, suffix, sep = ''), p, width=6, height=4)
  p
}

