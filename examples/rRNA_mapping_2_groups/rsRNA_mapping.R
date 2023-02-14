#!/user/bin/Rscript
library(reshape2)
library(ggplot2)
library(ggpubr)
library(wrMisc) # rowSEMs

experi=read.csv("experiments_design_for_deseq2Rscript.csv",header=TRUE)
Hsamples=experi[experi[,3] == "H",1]
Lsamples=experi[experi[,3] == "L",1]
Nsamples=experi[experi[,3] == "N",1]

#rsRNA mapping
data=NULL
for(i in list.files(".","28S")){  #list.files将文件夹中的文件名存到列表当中，可以用于批量导入文件
    s=gsub("_R1_001_28S_rRNA_mapping_revised.txt","",i)  #gsub("目标字符", "替换字符", 对象)
    d=read.table(i,header=TRUE,sep="\t")
    d$sample=s
    data=rbind(data,d)
    print(length(d$length))
}

d1=dcast(data,sample~length,value.var="RPM")

rownames(d1)=d1$sample

d1a=d1[,-1]

re.data=as.data.frame(t(d1a))

#specifc a specific region
region_start=135-20
region_end=135+19+20

re.data=re.data[as.character(region_start:region_end),]


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

p=ggplot(dp, aes(x = length, y = value, color=variable,group=variable))+geom_line(size=1)+theme_pubr()+
    geom_ribbon(aes(ymin = low, ymax = high,fill=variable), alpha = 0.1,colour = NA)+scale_color_manual(values = c("#00AFBB","#E7B800","#FC4E07"))+
    scale_fill_manual(values = c("#00AFBB","#E7B800","#FC4E07"))+ylab("Coverage(RPM)")+xlab("Length of nucleotides")+labs(color="group")+labs(fill="group")+
    scale_x_discrete(breaks = seq(min(dp$length), max(dp$length), by = 10))
p
# ggsave("RNA5_CGGGGCGCGGGACATGTGG.pdf",p,width=8,height=4)

