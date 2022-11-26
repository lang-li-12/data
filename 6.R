library(tidyverse) 
library(clusterProfiler)
library(org.Hs.eg.db)

library(readr)
MOCS2_DT_cor <- read_delim("GO/MOCS2-DT_cor.txt", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)

colnames(MOCS2_DT_cor) = 'Gene'

ego = enrichGO(
  gene = MOCS2_DT_cor$Gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego_res = ego@result
save(ego,ego_res,file = 'GO.Rdata')

barplot(ego,showCategory = 20,color = 'pvalue')

dotplot(ego,showCategory = 20)
