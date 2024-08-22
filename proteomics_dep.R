library(DEP) # 1.22.0
library(dplyr) # 1.1.4
setwd('.')
full_data <- read.delim('source_data/proteomics_MaxQuant_LFQ.txt',
                        stringsAsFactors = FALSE
)
exp_design <- read.delim('source_data/proteomics_dep_expdesign.txt',
                         stringsAsFactors = FALSE
)
exp_design <- exp_design[c(1:15),]
se <- import_MaxQuant(full_data, exp_design)
data_filt <- filter_missval(se, thr = 1)
data_norm <- normalize_vsn(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff_all_contrasts <- test_diff(data_imp, type = "all")
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))
data_results <- get_results(dep)
write.table(data_results,'source_data/dep_diffexp.txt',sep='\t')

nrow(data_results[(data_results$X40_lag_vs_X40_O_N_significant) & (data_results$X40_lag_vs_X40_O_N_ratio>0),])