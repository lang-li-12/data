#数据整理
library(tidyverse)
library(limma)

exp_693 = read.csv("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp_325 = read.csv("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp_lnc = inner_join(exp_325,exp_693,by = 'Gene_Name')
exp_lnc = column_to_rownames(exp_lnc,var = 'Gene_Name')

boxplot(exp_lnc,outline=FALSE,notch=T,las=2)  

sample_trait <- read_delim("sample_trait.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
exp_lnc = exp_lnc[,c(colnames(exp_lnc) %in% sample_trait$CGGA_ID)]
boxplot(exp_lnc,outline=FALSE,notch=T,las=2)


lncRNA <- read_delim("lncRNA.txt", delim = "\t", 
                     escape_double = FALSE, trim_ws = TRUE)
lncRNA = lncRNA[,-1,drop = F]
exp_lnc = rownames_to_column(exp_lnc,var = 'lncRNA')
exp_lnc = inner_join(exp_lnc,lncRNA,by = 'lncRNA')
exp_lnc = column_to_rownames(exp_lnc,var = 'lncRNA')
