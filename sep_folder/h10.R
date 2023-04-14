#作业十：
rm(list=ls())
options(stringsAsFactors = F)
#数据准备：下载数据集GSE42872的表达矩阵
f='GSE42872_eSet.Rdata'
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE42872_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
#获取分组列表
pd=pData(a)
#按照title分组
group_list=unlist(lapply(pd$title,function(x){strsplit(x," ")[[1]][4]}))

#进行差异分析---直接用代码
exprSet=dat
# DEG by limma
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design#样本分组，属于为1， 不属于为0
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels=design)
contrast.matrix#矩阵把progres组跟stable进行差异分析比较
##step1
fit <- lmFit(exprSet,design)
##step2
fit2<- contrasts.fit(fit, contrast.matrix) ##-#REE
fit2<- eBayes(fit2) ## default no trend !!!
##Bayes with trend-TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
#write.csv(nrDEGZ,"limma_notrend.results.csv", quote = F)
head(nrDEG)
dat=nrDEG

#画个火山图
library(ggpubr)
library(ggthemes)
nrDEG$logP=-log10(nrDEG$adj.P.Val)
head(nrDEG)
ggscatter(nrDEG,x="logFC",y="logP")+theme_base()
#火山图美化
#新加一列：
nrDEG$Group="not-significant"
nrDEG$Group[which((nrDEG$adj.P.Val<0.05) & (nrDEG$logFC>2))] <- "up_regulated"
nrDEG$Group[which((nrDEG$adj.P.Val<0.05) & (nrDEG$logFC< -2))] <- "down_regulated"
table(nrDEG$Group)#查看上调和下调基因数目
ggscatter(nrDEG,x="logFC",y="logP",color = "Group")+theme_base()
#改变火山图颜色(palette)和点的大小(size)
ggscatter(nrDEG,x="logFC",y="logP",color = "Group",
          palette = c("green","gray","red"),
          size=1)+theme_base()
#为火山图添加1ogP分界线(geom_hline)和1ogFC分界线(geom_vline)
library(ggplot2)
ggscatter(nrDEG,x="logFC",y="logP",color = "Group",
          palette = c("green","gray","red"),
          size=1)+theme_base()+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-2,2),linetype="dashed")
