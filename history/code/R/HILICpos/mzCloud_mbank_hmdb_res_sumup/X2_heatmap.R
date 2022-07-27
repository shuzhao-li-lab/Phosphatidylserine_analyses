# library
library(tidyr)
library(dplyr)
detach(package:plyr)

# read tables

mcg_res_fdr <- "../../output/RPneg/X4_mcg/1616888916.243331.RPneg_rm_2outliers_padj_default_m_negative_HEUvsHUU/"
lofCpd_df <- read.table(paste(mcg_res_fdr,"tables/ListOfEmpiricalCompounds.tsv",sep = ""), 
                        sep = "\t", header = 1, quote = "\"") # fill = TRUE, 

lofPathway_df <- read.table(paste(mcg_res_fdr,"tables/mcg_pathwayanalysis_RPneg_rm_2outliers_padj_default_m_negative_HEUvsHUU.tsv",sep = ""), 
                        sep = "\t", header = 1, quote = "\"") # fill = TRUE, 


featab <- read.csv("../../output/RPneg/X2.1_feat_tab_rm_outliers/rm_HUU_152_153_featureValues_summarized_log2scale.csv")

ttest_res_df <- read.table("../../output/RPneg/X3.1_ttest_0327_rm_HUU152_HUU153/ttest_res_full_report_HEUvsHUU.txt", header = 1)

output_dir <- "../../output/RPneg/X4_mcg/heatmaps/"
dir.create(output_dir)

threshold <- 0.2
adduct_vec <- c("M-H2O-H[-]","M-H[-]", "M-2H[2-]")

out_pdf_file = paste(output_dir ,"RPneg_Empd_based_on_stat_selected_adducts_","padj_",as.character(threshold),".pdf",sep = "")  #change the number
out_pdf_file

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
mm_df$EID <- factor(mm_df$EID)
mm_df$X <- factor(mm_df$X)

View(mm_df)
# https://stackoverflow.com/questions/26923862/why-are-my-dplyr-group-by-summarize-not-working-properly-name-collision-with
detach(package:plyr) # This is critical to avoid single-row summary which is used in plyr
EIDbased_sum_mm_df <- mm_df %>% group_by(EID) %>% summarize(pathway_collapsed = paste(pathway, collapse = ";"))
dim(mm_df)
dim(EIDbased_sum_mm_df)

path_vec <- c()
for (item in EIDbased_sum_mm_df$pathway_collapsed) {
  temp <- strsplit(item,";")[[1]]
  temp <- unique(temp)
  temp <- temp[!(temp %in% c("NA"))]
  if (length(temp) == 0) {
    temp <- ""
  } else {
    temp <- paste(temp, collapse = ";")
  }
  
  
  # temp <- gsub(";NA", "",item)
  # temp <- gsub("NA","",temp)
  path_vec <- c(path_vec,temp)
}

EIDbased_sum_mm_df$pathway_collapsed_clean <- path_vec

mmm_df <- merge(mm_df,  EIDbased_sum_mm_df, by = "EID")
mmm_df$pathway <- NULL
mmm_df <- unique(mmm_df)

maxLength = 50
# m_df$cbn_name <- substring(paste0(m_df$X, m_df$compound_names,sep = "_"),1,50)
mmm_df$cbn_name <- substring(paste(mmm_df$X, mmm_df$compound_names, sep = "|"),1,50)

duplicated_l <- mmm_df$cbn_name[duplicated(mmm_df$cbn_name)];duplicated_l
if(TRUE) {
  unique_mmm_df <- mmm_df[!duplicated(mmm_df$cbn_name),]
}


row.names(unique_mmm_df) <- unique_mmm_df$cbn_name
colnames(unique_mmm_df)

filt_m_df <- unique_mmm_df[unique_mmm_df$pval < threshold,] #
colnames(filt_m_df)
# [1] "EID"                     "X"                       "mzmed"                  
# [4] "mzmin"                   "mzmax"                   "rtmed"                  
# [7] "rtmin"                   "rtmax"                   "npeaks"                 
# [10] "HEU"                     "HUU"                     "ms_level"               
# [13] "HUU_082"                 "HUU_166"                 "HEU_234"                
# [16] "HEU_210"                 "HUU_169"                 "HEU_223"                
# [19] "HEU_188"                 "HEU_229"                 "HEU_190"                
# [22] "HUU_164"                 "HUU_161"                 "HUU_167"                
# [25] "HEU_078"                 "HUU_149"                 "HUU_081"                
# [28] "HUU_162"                 "HEU_075"                 "HEU_203"                
# [31] "HEU_232"                 "HUU_080"                 "HEU_231"                
# [34] "HUU_154"                 "HEU_195"                 "HEU_222"                
# [37] "HUU_168"                 "HEU_194"                 "HUU_155"                
# [40] "HEU_189"                 "HEU_191"                 "HUU_170"                
# [43] "HUU_160"                 "HEU_235"                 "HUU_158"                
# [46] "HUU_165"                 "HUU_171"                 "HEU_079"                
# [49] "HEU_067"                 "HEU_226"                 "t_score"                
# [52] "pval"                    "padj"                    "str_row_ion"            
# [55] "compounds"               "compound_names"          "adduct_type"            
# [58] "pathway_collapsed"       "pathway_collapsed_clean" "cbn_name"            

# Select only the one with good adduct type to present in the heatmap
filt_m_df <- filt_m_df[filt_m_df$adduct_type %in% adduct_vec,]

HT <- filt_m_df[,13:50] #
annot_row <- filt_m_df %>% select(adduct_type, pathway_collapsed_clean)
row.names(annot_row) <- row.names(filt_m_df)
colnames(annot_row) <- c("adduct", "pathways")


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
annotation_col = data.frame(label = sapply(colnames(HT),function(x)strsplit(x,"_")[[1]][1])) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(z_transformed) 

#annotation_row = data.frame(GeneClass = factor(HT$pathways)) # metabolism A...
#rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
label_color = c("#ff0000","#0000ff")
names(label_color) = c("HEU","HUU") #!
label_color
# HEU       HUU 
# "#ff0000" "#0000ff" 

range(z_transformed) #check the range of z_transformed
breaksList = seq(-2, 2, 0.1)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()


mat_breaks <- seq(min(z_transformed), max(z_transformed), length.out = 1)



pdf(file=out_pdf_file, width =100,height=30, paper = "special",onefile=FALSE)
par(mar = c(5,15,5,15))
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = annot_row, clustering_method = "complete", breaks = breaksList, 
         color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()