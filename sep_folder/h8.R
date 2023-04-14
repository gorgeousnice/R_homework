#作业八：
#找到 GPL6244 platform of Affymetrix Human Gene 1.0 ST Array 
#对应的R的bioconductor注释包:
options()$repos
options()$BioC_mirror 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#在这一步输入自己找到的R包
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)
options()$repos
options()$BioC_mirror
