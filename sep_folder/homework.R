#十题小作业汇总
rm(list=ls())
options(stringsAsFactors = F)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业一
dat <- read.table("ensembl.csv",head=F,sep="\t")
library(stringr)
dat$ensembl_id=str_split(dat$V1,"[.]",simplify = T)[,1]
library(org.Hs.eg.db)
g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)
bridge=merge(dat,g2e,by="ensembl_id")
new_data=merge(bridge,g2s,by="gene_id")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业二：
dat <- read.table("e2.txt",head=F,sep="\t")
colnames(dat)="probe_id"#改列名探针id
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
new_data=merge(dat,ids,by="probe_id")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业三：
#加载注释包，找到TP53所对应的基因
suppressPackageStartupMessages(library(hgu95av2.db))
ids = toTable(hgu95av2SYMBOL)
# 在ids的symbol列中搜索TP53
##方法一：可以直接View(ids)中搜索
##方法二： 
ids[ids$symbol %in% "TP53",]# %in% 返回逻辑值
ids[ids$symbol%in%"TP53",][,1] #返回TP53的探针名字

#获取表达矩阵与分组信息：
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
exprSet = exprs(sCLLex)# exprs()函数提取表达矩阵
# pData()函数获取分组信息:
pd = pData(sCLLex)# disease列:是按照疾病progressive进展-stable稳定分组
#绘制三个探针的boxplot
boxplot(exprSet["1939_at",] ~ pd$Disease) ##sig显著：有差异
boxplot(exprSet['1974_s_at',] ~ pd$Disease)
boxplot(exprSet['31618_at',] ~ pd$Disease)


#ggpubr绘制boxplot图形美化：
d <-cbind(as.data.frame(exprSet['1939_at',]),
          as.data.frame(pd$Disease))#合并1939_at的表达量和Disease
head(d)
suppressPackageStartupMessages(library(ggpubr))#加载包
colnames(d)<-c("Expression","Disease")#修改数据框列名
p<-ggboxplot(d, x="Disease", y ="Expression",
             palette = "jco",
             add = "jitter")
# 设置图形的参数
opar<-par(no.readonly=T)
# Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业四：
#在http://www.cbioportal.org/index.do网站
#找到BRCA1基因在TCGA数据库的乳腺癌数据集(Breast Invasive Carcinoma (TCGA, PanCancer Atlas))的表达情况
#Data type:clinical attribute---Subtype
dat <- read.table("plot_4.txt",head=T,fill=T,sep="\t")
head(dat)
colnames(dat)=c("Sample.Id ","Subtype","Expression","Mutations","CNA")
library(ggstatsplot)
ggbetweenstats(data=dat,x="Subtype",y="Expression")
library(ggplot2)
ggsave("Cbioportal_BRCA1_breast.png")#保存小提琴图
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业五：
#TP53基因在TCGA数据库的乳腺癌数据集
#按表达量分组绘制生存曲线http://www.oncolnc.org/
dat <- read.table("BRCA_7157_50_50.csv",head=T,fill=T,sep=",")
head(dat)
#重画网页中“生存分析”的图：
#以下代码是别人写好的，保证准备的数据中有Status,Expression,Days即可
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_TCGA.png")#保存为图片格式

#但发现p-value值没有小于0.05，并不显著：
#将作业四的数据plot_4.txt与BRCA_7157_50_50.csv合并分析
dat1 <- read.table("plot_4.txt",head=T,fill=T,sep="\t")
dat2 <- read.table("BRCA_7157_50_50.csv",head=T,fill=T,sep=",")       
colnames(dat1)=c("Patient","Subtype","Expression","Mutations","CNA")
dat1$Patient=substring(dat1$Patient,1,12)#12位以后的不要
merge_dat=merge(dat1,dat2,by="Patient")
table(merge_dat$Subtype)#查看Subtype类型
#分别按照类型分类，画生存分析图
#类型一：
dat=merge_dat[merge_dat$Subtype=="BRCA_Basal",]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_Basal_TCGA.png")#保存为图片格式
#类型二：
dat=merge_dat[merge_dat$Subtype=="BRCA_Her2",]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_Her2_TCGA.png")#保存为图片格式
#类型三：
dat=merge_dat[merge_dat$Subtype=="BRCA_LumA",]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_LumA_TCGA.png")#保存为图片格
#类型四：***显著
dat=merge_dat[merge_dat$Subtype=="BRCA_LumB",]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_LumB_TCGA.png")#保存为图片格
#类型五：
dat=merge_dat[merge_dat$Subtype=="BRCA_Normal",]
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status=="Dead",1,0)#拼写！
sfit <- survfit(Surv(Days,Status)~Group,data=dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)#无阴影，显示p-value
ggsave("survival_TP53_in_BRCA_Normal_TCGA.png")#保存为图片格式
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业六：
#获取数据
f='GSE17215_eSet.Rdata'
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE17215', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE17215_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)#获取数据集GSE17215的表达矩阵

#将GSE17215的表达矩阵中的探针id改成Symbol，缩小表达矩阵绘制热图
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
dat=dat[ids$probe_id,]#筛选删除没有对应Symbol的探针id的GSE17215的表达矩阵
dat[1:4,1:4] 
ids$median=apply(dat,1,median)#计算每一行的中位数
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]#去重复
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol#转换列名
dat[1:4,1:4]  
dim(dat)

#导入基因：
ng='ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T'
ng=strsplit(ng,' ')[[1]]
table(ng %in%  rownames(dat))
ng=ng[ng %in%  rownames(dat)]#过滤已经ng不在的dat里的
dat=dat[ng,]
dat=log2(dat)#数据标准化
pheatmap::pheatmap(dat,scale = 'row')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业七：
#获取数据：
f='GSE24673_eSet.Rdata'
library(GEOquery)

if(!file.exists(f)){
  gset <- getGEO('GSE24673', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE24673_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)

#获取分组数据
pd=pData(a)
group_list=c('rbc','rbc','rbc',
             'rbn','rbn','rbn',
             'rbc','rbc','rbc',
             'normal','normal')#创建根据source_name_ch1分组
dat[1:4,1:4]
M=cor(dat)
pheatmap::pheatmap(M)
tmp=data.frame(g=group_list)#加上分组信息
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M,annotation_col = tmp)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业八：
#找到 GPL6244 platform of Affymetrix Human Gene 1.0 ST Array 
                               #对应的R的bioconductor注释包:
options()$repos
options()$BioC_mirror 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)
options()$repos
options()$BioC_mirror
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业九：
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
#挑选出所有样本的(平均表达量/sd/mad/)最大的探针的列名
sort(apply(dat,1,mean),decreasing = T)[1]
sort(apply(dat,1,sd),decreasing = T)[1]
sort(apply(dat,1,mad),decreasing = T)[1]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业十：
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

pd=pData(a)
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
dim(dat)

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
