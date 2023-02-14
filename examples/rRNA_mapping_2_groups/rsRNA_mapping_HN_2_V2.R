#!/user/bin/Rscript
library(reshape2)
library(ggplot2)
library(ggpubr)
library(wrMisc) # rowSEMs


# experi=read.csv("experiments_design_for_deseq2Rscript.csv",header=TRUE)
# Hsamples=experi[experi[,3] == "T",1]
# Nsamples=experi[experi[,3] == "N",1]

rsRNAmapping <- function(treat_group) {
Hsamples=c(Hsamples)
Nsamples=c(Nsamples)
#rsRNA mapping
data=NULL
for(i in list.files(".","28S")){  #list.files将文件夹中的文件名存到列表当中，可以用于批量导入文件
  s=gsub("_28S_rRNA_mapping_revised.txt","",i)  #gsub("目标字符", "替换字符", 对象)
  d=read.table(i,header=TRUE,sep="\t")
  d$sample=s
  data=rbind(data,d)
  print(length(d$length))
}

d1=dcast(data,sample~length,value.var="RPM") #长数据变宽数据，～左边是行名，右边是列名，value.var：选择哪些列的数据转换

rownames(d1)=d1$sample

d1a=d1[,-1]

re.data=as.data.frame(t(d1a)) # 转至

#specifc a specific region
region_start=135-20  #重要变量
region_end=135+19+20  #重要变量
# region_start=1 #重要变量
# region_end=ncol(d1a)  #重要变量

re.data=re.data[as.character(region_start:region_end),]  #取符合的行


hsem=rowSEMs(re.data[,Hsamples]) # SEMs平均值的标准误差
# lsem=rowSEMs(re.data[,Lsamples])
nsem=rowSEMs(re.data[,Nsamples])

tmean=rowMeans(re.data[,Hsamples])  #Hsamples列名
# lmean=rowMeans(re.data[,Lsamples])
nmean=rowMeans(re.data[,Nsamples])

# re.data=data.frame(length=rownames(re.data),tmean=tmean,lmean=lmean,nmean=nmean)
re.data=data.frame(length=rownames(re.data),tmean=tmean,nmean=nmean)

dp=melt(re.data, id.vars = "length")  # 宽数据变长数据，id.vars：选择用来做主键的列

# dp$sem=c(hsem,lsem,nsem)
dp$sem=c(hsem,nsem)
dp$low=dp$value-dp$sem
dp$high=dp$value+dp$sem

dp$variable=factor(dp$variable,levels=c("tmean","nmean")) # "nmean",

p=ggplot(dp, aes(x = length, y = value, color=variable,group=variable))+
  geom_line(size=1)+theme_pubr()+
  geom_ribbon(aes(ymin = low, ymax = high,fill=variable), alpha = 0.1,colour = NA)+
  scale_color_manual(values = c("#00AFBB","#E7B800","#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB","#E7B800","#FC4E07"))+
  ylab("Coverage(RPM)")+xlab("Length of nucleotides")+labs(color="group")+labs(fill="group")+
  scale_x_discrete(breaks = seq(min(dp$length), max(dp$length), by = 10 ))  # by =10  (region_end-region_start+1)/6
p

# ggsave("RNA5_2_tn_rsRNA-249_CGGGGCGCGGGACATGTGG.pdf",p,width=8,height=4)
}