rm(list = ls())

library(plyr)

xcms_Featab <- read.csv("../../data/output/HILICpos_0402/cleanup_and_stat_test/featureValues_summarized_log2scale.csv")
head(xcms_Featab)

lf_val <- t(xcms_Featab[,2:ncol(xcms_Featab)]) # lf: long-formatted
colnames(lf_val) <- xcms_Featab$X
lf_val[1:5,1:10]
# FT0001   FT0002   FT0003   FT0004   FT0005   FT0006
# G2_R5pos_1 19.80967       NA 17.15374 17.21725 17.94974 22.93556
# G2_R5pos_2 19.80213       NA 24.08467 17.74124 17.96691 16.19684
# G2_R5pos_3 19.98046       NA 17.86672 17.79909 17.88582 23.29473
# G1_Naive_1 19.84599       NA 24.13974 17.97085 18.08051 17.83918
# G1_Naive_2 19.78193 18.55055 17.62319 18.24213 18.26202 19.20744

if(F) { # If log2 scale was not done
  lf_log2_val <- log(lf_val,2)
  lf_log2_val[1:5,1:10]
  lf_log2_val <- as.data.frame(lf_log2_val)
} else {
  lf_log2_val <- as.data.frame(lf_val)
}


#assignment of group/label
## make sure you have already turn mtx back to df. Otherwise, error may happen: Coercing LHS to a list
lf_log2_val$label <- factor(sapply(rownames(lf_log2_val), function(x)strsplit(x,"_")[[1]][2])) 

lf_log2_val$label
# You may want reorder the factor levels of the label before ttest.
lf_log2_val$label <- factor(lf_log2_val$label, levels = c("R5pos", "R5neg", "Naive"))

label_index = which(colnames(lf_log2_val)=="label") #or using grep

combn_ttest <- combn(unique(lf_log2_val$label),2)

list_data <- list()

for (i in 1:dim(combn_ttest)[2]) {
  temp_df <- lf_log2_val[lf_log2_val$label %in% combn_ttest[,i],]
  list_data[[i]]  <- lapply(temp_df[-label_index], function(x) {
    return(tryCatch(t.test(x ~ temp_df$label, var.equal = TRUE),error=function(e) NULL)) 
    # here I use default setting, no assumption of equal variance is made. And default is two.sided.
  })
}


#https://stackoverflow.com/questions/31468148/using-lapply-to-create-t-test-table/31468483
ttest_res_df_list = list() 
for (j in 1:length(list_data)) {
  
  data <- list_data[[j]]
  temp_row_list <- list()
  for (i in 1:length(data)) {
    if(is.null(data[i][[1]])) {
      temp_row_list[[i]] <- data.frame(names(data)[i],0,1)
      colnames(temp_row_list[[i]]) <- c("FT_ID","t_score","pval")
    } else {
      temp_row_list[[i]] <- data.frame(names(data)[i], data[[i]]$statistic, data[[i]]$p.value)
      colnames(temp_row_list[[i]]) <- c("FT_ID","t_score","pval")
    }
    
  }
  ttest_res_df = do.call(rbind, temp_row_list)
  ttest_res_df$padj <- p.adjust(ttest_res_df$pval, method = 'BH')
  ttest_res_df_list[[j]] <- ttest_res_df
}
 
#----------

# read the FT definition and prepare mummichog input; write one complete output & a mummichog input
featDef_df <- read.csv("../../data/output/HILICpos_0402/xcms/featureDefinitions.csv")
colnames(featDef_df)
# [1] "X"        "mzmed"    "mzmin"    "mzmax"    "rtmed"    "rtmin"    "rtmax"    "npeaks"   "G1_Naive" "G2_R5pos" "G3_R5neg"
# [12] "ms_level"

featVal_df <- xcms_Featab 
colnames(featVal_df)
#[1] "X"          "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2" "G3_R5neg_3"

featful_df <- merge(featDef_df, featVal_df, by = "X")
colnames(featful_df)
# [1] "X"          "mzmed"      "mzmin"      "mzmax"      "rtmed"      "rtmin"      "rtmax"      "npeaks"     "G1_Naive"   "G2_R5pos"  
# [11] "G3_R5neg"   "ms_level"   "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2"
# [21] "G3_R5neg_3"

colnames(ttest_res_df_list[[1]])
#[1] "FT_ID"   "t_score" "pval"    "padj"

m_featDef_ttest_df_list <- list()
for (i in 1:length(ttest_res_df_list)) {
  m_featDef_ttest_df_list[[i]] <- merge(featful_df, ttest_res_df_list[[i]], by.x = "X", by.y = "FT_ID")
}


## write the full report, for mummichog with either padj or raw pval

output_dir <- "../../data/output/HILICpos_0402/ttest_equal/"
dir.create(output_dir)

colnames(m_featDef_ttest_df_list[[1]])
# [1] "X"        "mzmed"    "mzmin"    "mzmax"    "rtmed"    "rtmin"    "rtmax"    "npeaks"   "G1_Naive" "G2_R5pos" "G3_R5neg"
# [12] "ms_level" "t_score"  "pval"     "padj"   

combn_ttest
# [,1]  [,2]  [,3] 
# [1,] R5pos R5pos Naive
# [2,] Naive R5neg R5neg
# Levels: R5pos R5neg Naive

for (i in 1:length(m_featDef_ttest_df_list)) {
  m_df <- m_featDef_ttest_df_list[[i]]
  
  A_str <- combn_ttest[,i][order(combn_ttest[,i])[1]]
  B_str <- combn_ttest[,i][order(combn_ttest[,i])[2]]
  AvsB_str <- paste(A_str,B_str, sep = "vs")
  
  write.table(m_df,paste(output_dir, "ttest_res_","full_report_", AvsB_str,".txt", sep = ""),
              sep = "\t",row.names = FALSE)
  
  essen_raw_pval_vec = c("mzmed","rtmed","pval","t_score","X","padj") # make sure columns are all in what you select
  essen_raw_pval_vec %in% colnames(m_df)
  
  full_col_name_raw_pval_vec <- c(essen_raw_pval_vec, colnames(m_df)[!(colnames(m_df) %in% essen_raw_pval_vec)])
  write.table(m_df[full_col_name_raw_pval_vec], 
              paste(output_dir, "ttest_res_","rawpval_",AvsB_str, ".txt", sep = ""), 
              sep = "\t",row.names = FALSE)
  
  essen_padj_vec = c("mzmed","rtmed","padj","t_score","X","pval") # make sure columns are all in what you select
  essen_padj_vec %in% colnames(m_df)
  
  full_col_name_padj_vec <- c(essen_padj_vec, colnames(m_df)[!(colnames(m_df) %in% essen_padj_vec)])
  
  write.table(m_df[full_col_name_padj_vec],
              paste(output_dir, "ttest_res_","padj_", AvsB_str,".txt", sep = ""), 
              sep = "\t",row.names = FALSE)

}
