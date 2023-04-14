#作业七：
rm(list=ls())
options(stringsAsFactors = F)
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
#计算样本的相关性并且绘制热图
dat[1:4,1:4]
M=cor(dat)
pheatmap::pheatmap(M)
tmp=data.frame(g=group_list)#加上分组信息
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M,annotation_col = tmp)