# remove the environment
rm(list = ls())

# library
library(tidyr)
library(dplyr)

# output directory create
output_dir <- "../../data/output/v0412_fix/HILICpos/heatmaps_based_on_mcg_res/"
dir.create(output_dir)

# Set mcg permutation p-value threshold
# mcg_pval_thres <- 0.05 # default
mcg_pval_thres <- 0.05

# Set the FT
FT_padj_threshold <- 0.05

only_abundant_adduct <- TRUE

key_name <- paste("HILICpos_R5posvsR5neg","mcgpadj",mcg_pval_thres ,"FTpadj",FT_padj_threshold,"adduct_abundant",only_abundant_adduct, sep = "_")


# read tables

mcg_res_fdr <- "../../data/output/v0412_fix/HILICpos/mcg_HILICpos/R5posvsR5neg/1618781765.861877.mcg_HILICpos_padj_0.05_R5posvsR5neg/"
lofCpd_df <- read.table(paste(mcg_res_fdr,"tables/ListOfEmpiricalCompounds.tsv",sep = ""), 
                        sep = "\t", header = 1, quote = "\"")

lofPathway_df <- read.table(paste(mcg_res_fdr,"tables/mcg_pathwayanalysis_mcg_HILICpos_padj_0.05_R5posvsR5neg.tsv",sep = ""), 
                        sep = "\t", header = 1, quote = "\"") 


featab <- read.csv("../../data/output/v0412_fix/HILICpos/feature_summary/featureValues_summarized_log2scale.csv")

ttest_res_df <- read.table("../../data/output/v0412_fix/HILICpos/ttest_equal/ttest_res_full_report_R5posvsR5neg.txt", header = 1)



# Examine the loaded tables

colnames(lofCpd_df)
# [1] "EID"              "massfeature_rows" "str_row_ion"      "compounds"        "compound_names"  

colnames(lofPathway_df)
# [1] "pathway"                         "overlap_size"                    "pathway_size"                   
# [4] "p.value"                         "overlap_EmpiricalCompounds..id." "overlap_features..id."          
# [7] "overlap_features..name." 

colnames(featab)
# [1] "X"          "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2" "G3_R5neg_3"

colnames(ttest_res_df)
# [1] "X"          "mzmed"      "mzmin"      "mzmax"      "rtmed"      "rtmin"      "rtmax"     
# [8] "npeaks"     "G1_Naive"   "G2_R5pos"   "G3_R5neg"   "ms_level"   "G2_R5pos_1" "G2_R5pos_2"
# [15] "G2_R5pos_3" "G1_Naive_1" "G1_Naive_2" "G1_Naive_3" "G3_R5neg_1" "G3_R5neg_2" "G3_R5neg_3"
# [22] "t_score"    "pval"       "padj" 

if(F) {  # If the ttest_res did not contain log2-scale vaules for inidivudal samples; do joining 
  ttest_res <- cbind(ttest_res_df,featab)
}

# Processing lofCpd_df
lofCpd_df <-tibble(lofCpd_df)
lofCpd_df_clean <- lofCpd_df %>% separate_rows(massfeature_rows, str_row_ion, sep = ";" , convert = TRUE)
lofCpd_df_clean$massfeature_rows <- gsub("row", "",lofCpd_df_clean$massfeature_rows)
lofCpd_df_clean$massfeature_rows <- formatC(
  sapply(lofCpd_df_clean$massfeature_rows,function(x)as.integer(x)), 
  width = 4, format = "d", flag = "0")
lofCpd_df_clean$massfeature_rows <- paste0("FT",lofCpd_df_clean$massfeature_rows)

lofCpd_df_clean$adduct_type <- sapply(lofCpd_df_clean$str_row_ion, function(x)strsplit(x,"_")[[1]][2])
lofCpd_df_clean[1:5,]
# EID   massfeature_rows str_row_ion   compounds          compound_names                      adduct_type
# <chr> <chr>            <chr>         <chr>              <chr>                               <chr>      
#   1 E4    FT0011           row11_M+H[1+] C00134             Putrescine                          M+H[1+]    
# 2 E6    FT0014           row14_M+H[1+] C00213;C00041;C00… Sarcosine; N-Methylglycine$L-Alani… M+H[1+]    
# 3 E7    FT0189           row189_M+Na[… C00213;C00041;C00… Sarcosine; N-Methylglycine$L-Alani… M+Na[1+]   
#                                       4 E15   FT0027           row27_M+H[1+] CE1950             cyanosulfurous acid anion           M+H[1+]    
#                                       5 E15   FT0030           row30_M+H[1+] CE1950             cyanosulfurous acid anion           M+H[1+] 

#Look at the lofPathway_df to determine the proper threshold and other pathways of interest.
# View(lofPathway_df) # You can drag to another monitor for better visualization


# Processing lofPathway_df
path_list_of_interest <- c("Glycolysis and Gluconeogenesis",
                           "Glycine, serine, alanine and threonine metabolism",
                           "Aspartate and asparagine metabolism",
                           "Arginine and Proline Metabolism",
                           "Glycerophospholipid metabolism",
                           "Methionine and cysteine metabolism",
                           "Tryptophan metabolism",
                           "Pentose phosphate pathway")
                           # "Leukotriene metabolism",
                           # "Fatty Acid Metabolism",
                           # "Glycerophospholipid metabolism",
                           # "Purine metabolism",
                           # "Pyrimidine metabolism")

lofPathway_df_clean <- lofPathway_df %>% separate_rows(overlap_EmpiricalCompounds..id., overlap_features..id., 
                                                       sep = "," , convert = TRUE) %>% 
  filter(p.value < mcg_pval_thres|pathway %in% path_list_of_interest) %>%
  select(overlap_EmpiricalCompounds..id., pathway) 
colnames(lofPathway_df_clean) <- c("EID","pathway")
lofPathway_df_clean[1:10,]


# Merge tables
m_df <- merge(ttest_res_df, lofCpd_df_clean, by.x = "X", by.y = "massfeature_rows")
mm_df <- merge(m_df,lofPathway_df_clean, by = "EID", all.x = TRUE)

mm_df$EID <- factor(mm_df$EID)

# Be careful of these stupid bug: 
# Don't load plyr second (after dplyr) or at all. The problem is that it's using plyr::summarise not dplyr::summarise
# Otherwise you will get only single row
# Be careful of this!

EIDbased_sum_mm_df <- mm_df %>% group_by(EID) %>% dplyr::summarise(pathway_collapsed = paste(pathway, collapse = ";"))

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
  
  
  temp <- gsub(";NA", "",item)
  temp <- gsub("NA","",temp)
  path_vec <- c(path_vec,temp)
}

EIDbased_sum_mm_df$pathway_collapsed_clean <- path_vec

mmm_df <- merge(mm_df, EIDbased_sum_mm_df, by = "EID")
mmm_df$pathway <- NULL
mmm_df$pathway_collapsed <- NULL
mmm_df <- unique(mmm_df)
dim(mmm_df)



# Set maximum length to avoid overload the page for both pathway and metabolite name
maxLength = 80
mmm_df$cbn_name <- substring(paste(mmm_df$X, mmm_df$compound_names,sep = "|"),1,maxLength)

## Add a choice to smplify the result.
if(FALSE) {
  mmm_df$pathway_collapsed_clean <- substring(mmm_df$pathway_collapsed_clean,1,maxLength)
} else {
  mmm_df$pathway_collapsed_clean <- sapply(mmm_df$pathway_collapsed_clean, function(x)strsplit(x, ";")[[1]][1])
}


row.names(mmm_df) <- mmm_df$cbn_name
colnames(mmm_df)

filt_m_df <- mmm_df[mmm_df$padj < FT_padj_threshold,]  # pval
colnames(filt_m_df)
# [1] "EID"                     "X"                       "mzmed"                  
# [4] "mzmin"                   "mzmax"                   "rtmed"                  
# [7] "rtmin"                   "rtmax"                   "npeaks"                 
# [10] "G1_Naive"                "G2_R5pos"                "G3_R5neg"               
# [13] "ms_level"                "G2_R5pos_1"              "G2_R5pos_2"             
# [16] "G2_R5pos_3"              "G1_Naive_1"              "G1_Naive_2"             
# [19] "G1_Naive_3"              "G3_R5neg_1"              "G3_R5neg_2"             
# [22] "G3_R5neg_3"              "t_score"                 "pval"                   
# [25] "padj"                    "str_row_ion"             "compounds"              
# [28] "compound_names"          "adduct_type"             "pathway_collapsed_clean"
# [31] "cbn_name"               

## Carefully check the samples you want to include in the heatmap. Each time will be different!

colnames(filt_m_df)[c(17:19,20:22)]
# [1] "G2_R5pos_1" "G2_R5pos_2" "G2_R5pos_3" "G3_R5neg_1" "G3_R5neg_2" "G3_R5neg_3"

HT <- filt_m_df[,c(17:19,20:22)] # Remeber to change accordingly!
annot_row <- filt_m_df %>% select(adduct_type, pathway_collapsed_clean)
row.names(annot_row) <- row.names(filt_m_df)
colnames(annot_row) <- c("adduct", "pathways")

adduct_vec <- c("M+H[1+]","M+Na[1+]")    # c("M-H2O-H[-]","M-H[-]", "M-2H[2-]")
annot_row$adduct[!(annot_row$adduct %in% adduct_vec)] = ""
if(only_abundant_adduct) {
  HT <- HT[(annot_row$adduct %in% adduct_vec),]
  annot_row <- annot_row[(annot_row$adduct %in% adduct_vec),]
}


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
annotation_col = data.frame(label = factor(c(rep("R5pos",3),rep("R5neg",3)))) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(z_transformed) 

#annotation_row = data.frame(GeneClass = factor(HT$pathways)) # metabolism A...
#rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
label_color = c("#0000ff","#ff0000")
names(label_color) = c("R5pos","R5neg") #!
label_color
# HEU       HUU 
# "#ff0000" "#0000ff" 

range(z_transformed) #check the range of z_transformed
breaksList = seq(-1.8, 1.8, 0.1)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()


mat_breaks <- seq(min(z_transformed), max(z_transformed), length.out = 1)

out_pdf_file = paste0(output_dir, key_name,".pdf")  #change the number


pdf(file=out_pdf_file, width =30,height=30, paper = "special",onefile=FALSE)
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = annot_row, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()