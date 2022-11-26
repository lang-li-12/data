#相关性分析
library(tidyverse)
library(limma)

exp_693 = read.csv("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp_325 = read.csv("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp = inner_join(exp_325,exp_693,by = 'Gene_Name')
exp = column_to_rownames(exp,var = 'Gene_Name')
library(readr)
sample_trait <- read_delim("sample_trait.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
exp = exp[,c(colnames(exp) %in% sample_trait$CGGA_ID)]
which(colnames(exp) == 'CGGA_703') #离群样本
exp = exp[,-42]

batch_cor <- function(gene){
  y <- as.numeric(exp[gene,])
  rownames <- rownames(exp)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exp[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

library(future.apply)
dd <- batch_cor("CTD-2366F13.1")
mRNA <- read_delim("mRNA.txt", delim = "\t",
                   escape_double = FALSE, trim_ws = TRUE)
colnames(mRNA) = 'mRNAs'
dd <- inner_join(dd,mRNA,by = 'mRNAs')

write.table(dd,file = "MOCS2-DT_cor.txt",sep = "\t",row.names = F,col.names = T,quote = F)
