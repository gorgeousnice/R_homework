#作业五：
rm(list=ls())
options(stringsAsFactors = F)
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

#但发现根据高低表达分组，p-value值没有小于0.05，并不显著：
#尝试按照Subtype细分组：
#将作业四的数据plot_4.txt与BRCA_7157_50_50.csv合并分析
dat1 <- read.table("plot_4.txt",head=T,fill=T,sep="\t")
dat2 <- read.table("BRCA_7157_50_50.csv",head=T,fill=T,sep=",")       
colnames(dat1)=c("Patient","Subtype","Expression","Mutations","CNA")
dat1$Patient=substring(dat1$Patient,1,12)#12位以后的不要
merge_dat=merge(dat1,dat2,by="Patient")

table(merge_dat$Subtype)#查看Subtype类型
#分别按照类型分类，画生存分析图
#类型一：Subtype类型为BRCA_Basal
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
ggsave("survival_TP53_in_BRCA_Normal_TCGA.png")#保存为图片格