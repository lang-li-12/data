library(readr)
exp <- read_delim("normalized_exp100736.txt",
                  delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE)
exp = exp[which(exp$...1 == 'CTD-2366F13.1'),]
exp = data.frame(MOCS2_DT = c( 6.77,7.08,7.14,8.44,8.01,8.40),
                 strain = c(1,1,1,0,0,0))

exp$strain = factor(exp$strain,levels = c(1,0),
                    labels = c('resistant','sensitive'))
levels(exp$strain)

colnames(exp)
var.test(MOCS2_DT ~ strain,data = exp)

t.test(MOCS2_DT ~ strain,data = exp,
       paired = F,var.equal = T)
