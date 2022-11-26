rm(list = ls())
dev.off()
#生存分析

library(readr)
deg_all <- read_delim("deg_all.txt", delim = "\t", 
                      escape_double = FALSE, trim_ws = TRUE)
pink <- read_delim("pink.txt", delim = "\t", 
                   escape_double = FALSE, trim_ws = TRUE)
library(tidyverse)
colnames(deg_all)[1] = 'lncRNA'
deg_lnc = inner_join(deg_all,pink,by = 'lncRNA')

logFC = 1   
P.Value = 0.05    

k1 = (deg_lnc$P.Value < P.Value)&(deg_lnc$logFC < -logFC)   
k2 = (deg_lnc$P.Value < P.Value)&(deg_lnc$logFC > logFC)     
deg_lnc$change = ifelse(k1,"down",ifelse(k2,"up","stable"))  

table(deg_lnc$change)
deg_lnc = deg_lnc[c(deg_lnc$change != 'stable'),]

sample_surv <- read_delim("sample_surv.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
exp_693 = read.csv("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp_325 = read.csv("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",sep = "\t",header=T,check.names=FALSE)
exp = inner_join(exp_325,exp_693,by = 'Gene_Name')
exp = column_to_rownames(exp,var = 'Gene_Name')
rm(exp_693)
rm(exp_325)

exp = as.data.frame(t(exp))
which(colnames(exp) %in% deg_lnc$lncRNA)
exp = exp[,which(colnames(exp) %in% deg_lnc$lncRNA)]
exp = rownames_to_column(exp,var = 'CGGA_ID')
exp_lncRNA = inner_join(exp,sample_surv,by = 'CGGA_ID')
exp_lncRNA = exp_lncRNA[,c(1,18,19,2:17)]
exp_lncRNA = column_to_rownames(exp_lncRNA,var = 'CGGA_ID')

which(exp_lncRNA$OS.time > 3600)
which(rownames(exp_lncRNA) %in% c('CGGA_703','CGGA_605','CGGA_1001'))
exp_lncRNA = exp_lncRNA[-c(1,23,60),]

#################################################################
library(survival)
library(survminer)

#使用平均数分组
m = mean(exp_lncRNA$`CTD-2366F13.1`)
exp_lncRNA$group = ifelse(exp_lncRNA$`CTD-2366F13.1` > m,"High","Low")
exp_lncRNA$group = factor(exp_lncRNA$group,levels = c("Low","High"))
class(exp_lncRNA$group)
table(exp_lncRNA$group)

#使用OS值
fitd = survdiff(Surv(OS.time,OS.status) ~ group,
                data = exp_lncRNA,
                na.action = na.exclude)
pValue = 1 - pchisq(fitd$chisq,length(fitd$n) - 1)

fit <- survfit(Surv(OS.time,OS.status) ~ group, data = exp_lncRNA)
ggsurvplot(fit, data = exp_lncRNA, risk.table = TRUE, conf.int = TRUE, pval = T)
