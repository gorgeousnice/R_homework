#作业六：
rm(list=ls())
options(stringsAsFactors = F)
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
ng=ng[ng %in%  rownames(dat)]#过滤dat里没有的ng
dat=dat[ng,]
dat=log2(dat)#数据标准化
pheatmap::pheatmap(dat,scale = 'row')
