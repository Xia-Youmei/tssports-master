#!/user/bin/Rscript
#tRNA mapping plot for comparsion group
#self.R experimentfile datafile
#' prepare_tsrna_mapping_plotdata
#'
#' @param tsRNAname  specific tsRNA name, which you can find in "_tRNA_mature_mapping.txt".
#' @param pattern  A uniform suffix for the filename results returned by the SPORTS1.1 tool.
#'
#' @return combined specific tsRNA mapping results, form 1 new file "specific tsRNA"_mappingplot.txt.
#' @export
#'
#' @examples tRNA_mapping_comparision(filename = "Homo_sapiens_tRNA-Glu-TTC-1_mappingplot.txt",
#' Hsamples=c("XWH5_S7_L004", "XWH4_S6_L004", "XWH1_S0_L000", "XWH2_S6_L003"),
#' Lsamples=c("XWH3_S5_L004",  "XWL4_S11_L004", "XWL1_S8_L004",  "XWL3_S0_L000",  "XWL2_S9_L004",  "XWL5_S1_L001" ),
#' Nsamples=c("XWN1_S1_L004", "XWN2_S2_L004", "XWN3_S5_L003")
#' )
#'
tRNA_mapping_comparision<- function(filename,
                                    Tsamples= FALSE,
                                    Hsamples= FALSE,
                                    Lsamples= FALSE,
                                    Nsamples,
                                    pattern = "_tRNA_mature_mapping.txt"){
  library(wrMisc)
  library(reshape2)
  library(ggpubr)

  if( Tsamples != FALSE ){
    Tsamples=Tsamples
    Nsamples=Nsamples

    tRNA.mapping <- read.table(filename ,sep="\t")

    tRNA.number <- nrow(tRNA.mapping)

    re.data=data.frame(length=1:length(strsplit(as.character(tRNA.mapping[1,2]), ",")[[1]]))
    for (i in 1:tRNA.number) {
      rpm <- strsplit(as.character(tRNA.mapping[i,2]), ",")
      title.name <- as.character(tRNA.mapping[i,1])
      len <- length(rpm[[1]])
      length = c(1:len)
      # print(len)
      length.data <- data.frame(RPM = as.numeric(rpm[[1]]))
      colnames(length.data)[1]=tRNA.mapping[i,3]
      re.data=cbind(re.data,length.data)
    }

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

    tRNA.mapping <- read.table(filename ,sep="\t")

    tRNA.number <- nrow(tRNA.mapping)

    re.data=data.frame(length=1:length(strsplit(as.character(tRNA.mapping[1,2]), ",")[[1]]))
    for (i in 1:tRNA.number) {
      rpm <- strsplit(as.character(tRNA.mapping[i,2]), ",")
      title.name <- as.character(tRNA.mapping[i,1])
      len <- length(rpm[[1]])
      length = c(1:len)
      # print(len)
      length.data <- data.frame(RPM = as.numeric(rpm[[1]]))
      colnames(length.data)[1]=tRNA.mapping[i,3]
      re.data=cbind(re.data,length.data)
    }

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
    scale_color_manual(values = c("#00AFBB","#FC4E07","#E7B800"))+
    scale_fill_manual(values = c("#00AFBB","#FC4E07","#E7B800"))+
    ylab("Coverage(RPM)")+xlab("Length of nucleotides")+labs(color="group")+labs(fill="group")
  p

  ggsave("Homo_sapiens_tRNA-Glu-TTC-T-N.pdf",p,width=8,height=4)
}

