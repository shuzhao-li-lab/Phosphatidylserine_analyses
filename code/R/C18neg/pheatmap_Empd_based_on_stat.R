# library
library(tidyr)
library(dplyr)

# read tables

mcg_res_fdr <- "../../data/output/RPneg_031821run/mcg/mcg_C18neg_equal_var_ttest/R5posvsR5neg/1616377073.842988.mcg_C18neg_padj_default_R5posvsR5neg/"
lofCpd_df <- read.table(paste(mcg_res_fdr,"tables/ListOfEmpiricalCompounds_mod.txt",sep = ""), 
                        sep = "\t", header = 1) # fill = TRUE, 

lofPathway_df <- read.table(paste(mcg_res_fdr,"tables/mcg_pathwayanalysis_mcg_C18neg_padj_default_R5posvsR5neg.tsv",sep = ""), 
                        sep = "\t", header = 1) # fill = TRUE, 


featab <- read.csv("../../data/output/RPneg_031821run/cleanup_and_stat_test/featureValues_summarized_log2scale.csv")

ttest_res_df <- read.table("../../data/output/RPneg_031821run/cleanup_and_stat_test/ttest_equal_variance/ttest_res_full_report_R5posvsR5neg.txt", header = 1)


colnames(lofCpd_df)
# [1] "EID"              "massfeature_rows" "str_row_ion"      "compounds"        "compound_names"  

colnames(lofPathway_df)
# [1] "pathway"                         "overlap_size"                    "pathway_size"                   
# [4] "p.value"                         "overlap_EmpiricalCompounds..id." "overlap_features..id."          
# [7] "overlap_features..name." 

colnames(featab)
# [1] "X"          "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2" "G3_R5neg_3"

colnames(ttest_res_df)
# [1] "X"          "mzmed"      "mzmin"      "mzmax"      "rtmed"      "rtmin"      "rtmax"      "npeaks"     "G1_Naive"   "G2_R5pos"  
# [11] "G3_R5neg"   "ms_level"   "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2"
# [21] "G3_R5neg_3" "t_score"    "pval"       "padj"    

if(F) {  # If the ttest_res did not contain log2-scale vaules for inidivudal samples; do joining 
  ttest_res <- cbind(ttest_res_df,featab)
}

#processing lofCpd_df
lofCpd_df <-tibble(lofCpd_df)
lofCpd_df_clean <- lofCpd_df %>% separate_rows(massfeature_rows, str_row_ion, sep = ";" , convert = TRUE)
lofCpd_df_clean$massfeature_rows <- gsub("row", "FT",lofCpd_df_clean$massfeature_rows)
lofCpd_df_clean$adduct_type <- sapply(lofCpd_df_clean$str_row_ion, function(x)strsplit(x,"_")[[1]][2])

#processing lofPathway_df
lofPathway_df_clean <- lofPathway_df %>% separate_rows(overlap_EmpiricalCompounds..id., overlap_features..id., 
                                                       sep = "," , convert = TRUE) %>% 
  filter(p.value < 0.05) %>%
  select(overlap_EmpiricalCompounds..id., pathway) 
colnames(lofPathway_df_clean) <- c("EID","pathway")

m_df <- merge(ttest_res_df, lofCpd_df_clean, by.x = "X", by.y = "massfeature_rows")
mm_df <- merge(m_df,lofPathway_df_clean, by = "EID", all.x = TRUE)
test <- mm_df %>%
          group_by(EID) %>%
          summarize(pathway_collapsed = paste(pathway, collapse = ";"))

path_vec <- c()
for (item in test$pathway_collapsed) {
  temp <- gsub(";NA", "",item)
  temp <- gsub("NA","",temp)
  path_vec <- c(path_vec,temp)
}
test$pathway_collapsed_clean <- path_vec

mmm_df <- merge(mm_df, test, by = "EID")
mmm_df$pathway <- NULL
mmm_df <- unique(mmm_df)


maxLength = 50
# m_df$cbn_name <- substring(paste0(m_df$X, m_df$compound_names,sep = "_"),1,50)
mmm_df$cbn_name <- paste(mmm_df$X, mmm_df$compound_names, sep = "|")

row.names(mmm_df) <- mmm_df$cbn_name
colnames(mmm_df)

filt_m_df <- mmm_df[mmm_df$pval < 0.05,] #
colnames(filt_m_df)
# [1] "EID"                     "X"                       "mzmed"                   "mzmin"                  
# [5] "mzmax"                   "rtmed"                   "rtmin"                   "rtmax"                  
# [9] "npeaks"                  "G1_Naive"                "G2_R5pos"                "G3_R5neg"               
# [13] "ms_level"                "G1_Naive_1"              "G1_Naive_2"              "G1_Naive_3"             
# [17] "G2_R5pos_1"              "G2_R5pos_2"              "G2_R5pos_3"              "G3_R5neg_1"             
# [21] "G3_R5neg_2"              "G3_R5neg_3"              "t_score"                 "pval"                   
# [25] "padj"                    "str_row_ion"             "compounds"               "compound_names"         
# [29] "X.y"                     "adduct_type"             "pathway_collapsed"       "pathway_collapsed_clean"
# [33] "cbn_name"            

HT <- filt_m_df[,17:22] #
annot_row <- filt_m_df %>% select(adduct_type, pathway_collapsed_clean)
row.names(annot_row) <- row.names(filt_m_df)
colnames(annot_row) <- c("adduct", "pathways")

adduct_vec <- c("M-H2O-H[-]","M-H[-]", "M-2H[2-]")
annot_row$adduct[!(annot_row$adduct %in% adduct_vec)] = ""

log2HT <- log(HT,2) #

#####
library(pheatmap)
library(RColorBrewer)

zscore <- function(x) {
  z <- (x - rowMeans(x))/ apply(x,1,sd)
  return(z)
}


#numeric_input <- log2HT[,1:ncol(log2HT)]  #!
z_transformed <-zscore(log2HT)  #log(numeric_input+1,2)
z_transformed <- z_transformed[complete.cases(z_transformed), ]

# Add annotation as described above, and change the name of annotation
annotation_col = data.frame(label = factor(c(rep("R5_pos",3),rep("R5neg",3)))) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(z_transformed) 

#annotation_row = data.frame(GeneClass = factor(HT$pathways)) # metabolism A...
#rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
label_color = c("#0000ff","#ff0000")
names(label_color) = c("R5_pos","R5neg") #!
label_color
# HEU       HUU 
# "#ff0000" "#0000ff" 

range(z_transformed) #check the range of z_transformed
breaksList = seq(-1.8, 1.8, 0.1)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()


mat_breaks <- seq(min(z_transformed), max(z_transformed), length.out = 1)

out_pdf_file = paste("../../data/output/RPneg_031821run/heatmaps_based_on_mcg_res/R5posvsR5neg_pval05",".pdf",sep = "")  #change the number


pdf(file=out_pdf_file, width =30,height=10, paper = "special",onefile=FALSE)
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = annot_row, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()