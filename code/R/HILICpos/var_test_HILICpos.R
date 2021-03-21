# Variance check

library(plyr)

xcms_Featab <- read.csv("../../../Rafi_Ahmed_exhaustedCD8Tcell/data/output/HILICpos_031821run/cleanup_and_stat_test/featureValues_summarized_log2scale.csv")
head(xcms_Featab)

lf_val <- t(xcms_Featab[,2:ncol(xcms_Featab)]) # lf: long-formatted
colnames(lf_val) <- xcms_Featab$X
lf_val[1:5,1:10]
# FT0001   FT0002   FT0003   FT0004   FT0005   FT0006   FT0007   FT0008   FT0009   FT0010
# G1_Naive_1 16.62220 17.55000 20.79110 19.71521 16.32105 20.16929       NA 23.28054 19.21986 14.50632
# G1_Naive_2 16.15796 17.59921 20.57279 20.22519 15.95561 20.64719 16.86536 23.90663       NA 14.57523
# G1_Naive_3 16.38336 17.78636 20.92688 20.25349 16.13105 20.67444       NA 23.65352 18.74238 14.49966
# G2_R5pos_1 17.27595 17.89888 21.28560 18.90236 17.07683 19.91885       NA 22.78709 18.53189 15.01651
# G2_R5pos_2 16.92078 17.89240 20.79049 19.38063 16.64326 19.84517       NA 22.98639       NA       NA

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

combn_test <- combn(unique(lf_log2_val$label),2)

for (i in 1:dim(combn_test)[2]) {
  temp_df <- lf_log2_val[lf_log2_val$label %in% combn_test[,i],]
  list_data[[i]]  <- lapply(temp_df[-label_index], function(x) {
    return(tryCatch(var.test(x ~ temp_df$label),error=function(e) NULL))
  })
}


#https://stackoverflow.com/questions/31468148/using-lapply-to-create-t-test-table/31468483
var_res_df_list = list() 
for (j in 1:length(list_data)) {
  
  data <- list_data[[j]]
  temp_row_list <- list()
  for (i in 1:length(data)) {
    if(is.null(data[i][[1]])) {
      temp_row_list[[i]] <- data.frame(names(data)[i],0,1)
      colnames(temp_row_list[[i]]) <- c("FT_ID","ratio","pval")
    } else {
      temp_row_list[[i]] <- data.frame(names(data)[i], data[[i]]$statistic, data[[i]]$p.value)
      colnames(temp_row_list[[i]]) <- c("FT_ID","ratio","pval")
    }
    
  }
  var_res_df = do.call(rbind, temp_row_list)
  var_res_df$padj <- p.adjust(var_res_df$pval, method = 'BH')
  var_res_df_list[[j]] <- var_res_df
}

hist(var_res_df_list[[1]]$pval)
hist(var_res_df_list[[1]]$pval)

# read the FT definition and prepare mummichog input; write one complete output & a mummichog input
featDef_df <- read.csv("../../data/output/HILICpos_031821run/featureDefinitions.csv")
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

colnames(var_res_df_list[[1]])
#[1] "FT_ID"   "ratio" "pval"    "padj"

m_featDef_vartest_df_list <- list()
for (i in 1:length(var_res_df_list)) {
  m_featDef_vartest_df_list[[i]] <- merge(featful_df, var_res_df_list[[i]], by.x = "X", by.y = "FT_ID")
}




## write the full report, for mummichog with either padj or raw pval
output_dir <- "../../data/output/HILICpos_031821run/cleanup_and_stat_test/"

combn_test
# [,1]  [,2]  [,3] 
# [1,] Naive Naive R5pos
# [2,] R5pos R5neg R5neg
# Levels: Naive R5neg R5pos

for (i in 1:length(m_featDef_vartest_df_list)) {
  m_df <- m_featDef_vartest_df_list[[i]]
  
  A_str <- combn_ttest[,i][order(combn_ttest[,i])[1]]
  B_str <- combn_ttest[,i][order(combn_ttest[,i])[2]]
  AvsB_str <- paste(A_str,B_str, sep = "vs")
  
  write.table(m_df,paste(output_dir, "var_test_res_","full_report_", AvsB_str,".txt", sep = ""),
              sep = "\t",row.names = FALSE)
  pdf(paste(output_dir, "var_test_res_","hist_", AvsB_str,".pdf", sep = ""))
  hist(var_res_df_list[[i]]$pval, breaks = 20)
  hist(var_res_df_list[[i]]$ratio, xlim = range(0,4),breaks = 500000)
  dev.off()
}

