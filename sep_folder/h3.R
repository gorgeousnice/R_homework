#作业三：
rm(list=ls())
options(stringsAsFactors = F)
#加载注释包，找到TP53所对应的探针名
suppressPackageStartupMessages(library(hgu95av2.db))
ids = toTable(hgu95av2SYMBOL)
##方法一：可以直接View(ids)中搜索
##方法二： 
ids[ids$symbol %in% "TP53",]# %in% 返回逻辑值
ids[ids$symbol%in%"TP53",][,1] #返回TP53的探针名字

#加载R包CLL内置的数据集的表达矩阵
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
          as.data.frame(pd$Disease))#合并1939_at各样本表达量信息和Disease分组
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

