#作业九：
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

#挑选出所有样本的(平均表达量/sd/mad/)最大的探针的列名
sort(apply(dat,1,mean),decreasing = T)[1]
sort(apply(dat,1,sd),decreasing = T)[1]
sort(apply(dat,1,mad),decreasing = T)[1]
