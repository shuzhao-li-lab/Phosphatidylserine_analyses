#Correct the xdata 

rm(list = ls())

library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)

output_dir = "../../data/output/HILICpos_0331analysis/xdata_corrected/"
dir.create(output_dir)

xdata <- readRDS(file = "../../data/output/HILICpos_031821run/xcms_old_run/HILICpos_xcms.rds")

basename(fileNames(xdata))
# [1] "G1_Naive_1.mzML"  "G1_Naive_2.mzML"  "G1_Naive_3.mzML"  "G2_R5pos_1.mzML" 
# [5] "G2_R5pos_2.mzML"  "G2_R5pos_3.mzML"  "G3_R5neg_1a.mzML" "G3_R5neg_1b.mzML"
# [9] "G3_R5neg_2a.mzML" "G3_R5neg_2b.mzML" "G3_R5neg_3a.mzML" "G3_R5neg_3b.mzML"

xdata$sample_group
# [1] "G2_R5pos" "G2_R5pos" "G2_R5pos" "G1_Naive" "G1_Naive" "G1_Naive" "G3_R5neg"
# [8] "G3_R5neg" "G3_R5neg" "G3_R5neg" "G3_R5neg" "G3_R5neg"

xdata$sample_name
# [1] "G2_R5pos_1"  "G2_R5pos_2"  "G2_R5pos_3"  "G1_Naive_1"  "G1_Naive_2" 
# [6] "G1_Naive_3"  "G3_R5neg_1a" "G3_R5neg_1b" "G3_R5neg_2a" "G3_R5neg_2b"
# [11] "G3_R5neg_3a" "G3_R5neg_3b"

# As can be seen here; xdata's phenodata doesn't match with basename
# To avoid any chromatogram extraction problem reshuffle the xdata's phenodata.


xdata$sample_group <- c(rep("G1_Naive",3), rep("G2_R5pos",3), rep("G3_R5neg",6))
xdata$sample_name <- c("G1_Naive_1","G1_Naive_2","G1_Naive_3",
                       "G2_R5pos_1","G2_R5pos_2","G2_R5pos_3",
                       "G3_R5neg_1a", "G3_R5neg_1b", 
                       "G3_R5neg_2a", "G3_R5neg_2b",
                       "G3_R5neg_3a", "G3_R5neg_3b")

saveRDS(xdata,paste0(output_dir, "pd_corr_HILICpos_xcms.rds"))


## Extract the features and log2 transform them
res <- quantify(xdata, value = "into")
ft_ints <- log2(assay(res, "raw")) 

#assign group_colors
group_colors <- c("#999999","#0000ff","#ff0000")
names(group_colors) <- unique(xdata$sample_group) 
group_colors
# G1_Naive  G2_R5pos  G3_R5neg 
# "#999999" "#0000ff" "#ff0000" 

## Perform the PCA omitting all features with an NA in any of the
## samples. Also, the intensities are mean centered.
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)

## Plot the PCA
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2],  col = "darkgrey", labels = "", #labels = xdata$sample_name
     pos = 3, cex = 2)

# print the pdf
pdf(paste(output_dir,"PCA_with_labels.pdf", sep = ""))
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2],  col = "darkgrey", labels = xdata$sample_name, #labels = "",
     pos = 1, cex = 1)
dev.off()