rm(list=ls())
options(stringsAsFactors = F)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#作业一:
dat <- read.table("ensembl.csv",head=F,sep="\t")
library(stringr)
dat$ensembl_id=str_split(dat$V1,"[.]",simplify = T)[,1]#点分割
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