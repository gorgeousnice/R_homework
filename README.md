# R_homework
来自曾老师的“R语言小作业-中级http://www.bio-info-trainee.com/3750.html”

可能你也跟我同样好奇，这些分析可以解决什么样的生物学问题；

其中homework_R是作业汇总，还有分开的代码在sep_folder文件夹中；

## 各种数据库简介： 
### GEO数据库(Gene Expression Ominibus)：基因表达综合数据库；  
包含肿瘤和非肿瘤；收录并整理了全球范国内研究工作者上传的微阵列芯片、二代测序以及其他形式的高通量基因组数据；
### TCGA数据库(The Cancer Genome Atlas)：癌症基因组图谱数据库
该数据库收录了33种癌症的20000多个样本的多种数据，包括了转录组表达数据、基因组变异数据、甲基化数据、临床数据等。
### cbioportal数据库：在cancer discovery上的对TCGA等数据库的基因组（包括DNA和RNA）分析工具

## 作业代码对应生物学问题:
### h1&2:基因名之间的转换：https://www.jianshu.com/p/3b27c32fa392
Ensembl_ID：Ensembl数据库中对基因的命名，以ENS开头； 
Gene_ID :是NCBI使用的能够对众多数据库进行联合搜索的搜索引擎, 其对不同的Gene进行了编号,就是一串数字;  
Gene_Symbol：是HGNC数据库为基因提供的官方名称，主要是按基因的功能起的名字，字母一般为英文全称的缩写，由大写字母和数字组成；  
probeset_id:探针id  
tips：这里有一篇关于GEO数据库ID转化的文章：https://www.bilibili.com/read/cv14560979

### h3:根据CLL包，绘制TP53基因表达的箱线图和美化
CLL包：包含慢性淋巴细胞白血病(CLL)基因表达数据。CLL数据中有24个样本，根据疾病进展情况被分类为进展或稳定progress-stable。 
根据progress-stable分组画箱线图：看该基因的表达量是否与疾病进展情况有关？？？

### h4:BRCA1基因在TCGA数据库的乳腺癌数据集的表达情况，画小提琴图：
从http://www.cbioportal.org/index.do进入：选择癌症类型和数据集

### h5:TP53基因在TCGA数据库的乳腺癌数据集的表达量分组的生存分析：
#### 关于生存分析曲线：
Time:生存时间or最后一次能够联系到这个病人的时间；随着时间推移，就会有病人该病去世    
Survival probability:每去世一个病人，整体的生存概率就会下降   
p-value：看两组病人的生存时间是否存在显著差异     
短线：最后一次能联系到病人的时间；病人不是由于癌症造成的死亡，联系不到了  

### h6:在数据集GSE17215的表达矩阵提取基因画热图 
包含将GSE17215的表达矩阵中的探针id改成Symbol***
这个热图可以得到什么：

### h7:在数据集GSE24673的表达矩阵计算样本的相关性并且绘制热图，并标记上样本分组信息
dat=exprs(a)获取样本数据 
pd=pData(a)获取分组数据 
这个相关性热图可以得到什么：

### h8:找到对应的R的bioconductor注释包 
Bioconductor开发的物种注释包系列集合了一个物种不同来源的注释信息，能够根据基因ID对其进行多种来源的注释，比如说基因的别名，基因的功能等。  
http://www.bio-info-trainee.com/1399.html

### h9:在数据集GSE42872的表达矩阵，挑选出所有样本的(平均表达量/sd/mad/)最大的探针

### h10：在数据集GSE42872的表达矩阵，根据分组使用limma做差异分析，得到差异结果矩阵画火山图













