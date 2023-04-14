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
