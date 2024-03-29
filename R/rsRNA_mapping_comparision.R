#' @title Plot rsRNA mapping comparision coverage profile
#'
#' @description
#'  Plot the coverage profile of the rsRNA, Solid curves presenting average RPM of 2 groups(treat VS control) or 3 groups(high VS control; low VS control). RPM, reads per million; nt, nucleotide. Colored band stood for 95% CI.
#'
#' @param Tsamples A collection of treat group file name prefixes.
#' @param Hsamples A collection of high group file name prefixes.
#' @param Lsamples A collection of low group file name prefixes.
#' @param Nsamples A collection of control group file name prefixes.
#' @param colors The color of the line in the line chart.
#' @param suffix A format suffix for a graph file, it can be jpg/png/tif/gif and so on.
#' @param rsRNAName Specific rsRNA name, like 28S/5S/5.8S and so on.
#' @param rsRNAFilePrefixJunctionSymbol RsRNA file prefix junction symbol, it is the concatenator that removes the Tsamples/Hsamples/Lsamples/Nsamples prefix and suffix.
#' @param rsRNAFileSuffix RsRNA file suffix, it's usually after rsRNAName.
#' @param region_start Length of nucleotides start region site, end region always add 59 based on start region site.
#'
#' @return the coverage profile of the rsRNA and form a new file "specific rsRNA".pdf.
#' @export
#'
#' @examples
#' 2 group : treat VS control
#' rsRNAmapping(Tsamples=c("XWH5_S7_L004", "XWH4_S6_L004", "XWH1_S0_L000", "XWH2_S6_L003"),
#' Nsamples=c("XWN1_S1_L004", "XWN2_S2_L004", "XWN3_S5_L003")
#'
#' 3 group : high VS control; low VS control
#' rsRNAmapping(Hsamples=c("XWH5_S7_L004", "XWH4_S6_L004", "XWH1_S0_L000", "XWH2_S6_L003"),
#' Lsamples=c("XWH3_S5_L004",  "XWL4_S11_L004", "XWL1_S8_L004",  "XWL3_S0_L000",  "XWL2_S9_L004",  "XWL5_S1_L001" ),
#' Nsamples=c("XWN1_S1_L004", "XWN2_S2_L004", "XWN3_S5_L003"),
#' rsRNAFilePrefixJunctionSymbol = '_R1_001_')
#'
#'
#'
rsRNAmapping<- function(Tsamples= NULL,
                        Hsamples= NULL,
                        Lsamples= NULL,
                        Nsamples,
                        rsRNAFilePrefixJunctionSymbol = '_',
                        rsRNAName = "28S",
                        rsRNAFileSuffix = '_rRNA_mapping_revised.txt',
                        region_start = 115,
                        colors = c("#00AFBB","#FC4E07","#E7B800"),
                        suffix='.pdf'){

  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(wrMisc)

  data=NULL
  for(i in list.files(".",rsRNAName)){
    rsRNANamePaste <- paste(rsRNAFilePrefixJunctionSymbol, rsRNAName, rsRNAFileSuffix, sep = '')
    s=gsub(rsRNANamePaste,"",i)
    d=read.table(i,header=TRUE,sep="\t")
    d$sample=s
    data=rbind(data,d)
    # print(length(d$length))
  }

  d1=dcast(data,sample~length,value.var="RPM")

  rownames(d1)=d1$sample

  d1a=d1[,-1]

  re.data=as.data.frame(t(d1a))

  #specifc a specific region
  region_start=region_start  #重要变量
  region_end=region_start+59  #重要变量

  re.data=re.data[as.character(region_start:region_end),]


  if( length(Tsamples) != 0 ){
    Tsamples=Tsamples
    Nsamples=Nsamples

    tsem=rowSEMs(re.data[,Tsamples])
    nsem=rowSEMs(re.data[,Nsamples])

    tmean=rowMeans(re.data[,Tsamples])
    nmean=rowMeans(re.data[,Nsamples])

    re.data=data.frame(length=rownames(re.data),tmean=tmean,nmean=nmean)

    dp=melt(re.data, id.vars = "length")  # 宽数据变长数据，id.vars：选择用来做主键的列

    dp$sem=c(tsem,nsem)
    dp$low=dp$value-dp$sem
    dp$high=dp$value+dp$sem

    dp$variable=factor(dp$variable,levels=c("tmean","nmean")) # "nmean",
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

    re.data=data.frame(length=rownames(re.data),hmean=hmean,lmean=lmean,nmean=nmean)

    dp=melt(re.data, id.vars = "length")

    dp$sem=c(hsem,lsem,nsem)
    dp$low=dp$value-dp$sem
    dp$high=dp$value+dp$sem

    dp$variable=factor(dp$variable,levels=c("nmean","hmean","lmean"))

  }
  p=ggplot(dp, aes(x = length, y = value, color=variable,group=variable))+
    geom_line(size=1)+theme_pubr()+
    geom_ribbon(aes(ymin = low, ymax = high,fill=variable), alpha = 0.1,colour = NA)+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    ylab("Coverage(RPM)")+xlab("Length of nucleotides")+labs(color="group")+labs(fill="group")+
    scale_x_discrete(breaks = seq(min(dp$length), max(dp$length), by = 10 ))  # by =10  (region_end-region_start+1)/6
  ggsave(paste(rsRNAName, '_rRNA_mapping',suffix, sep = ''), p, width=6, height=4)
  p
}
