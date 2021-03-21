# This is to collapse the technical replicates occur in some samples.

library(corrplot)

output_dir = "../../data/output/HILICpos_031821run/cleanup_and_stat_test/"
if(T) {
  dir.create(output_dir)
}

df <- read.csv("../../data/output/HILICpos_031821run/featureValues.csv", row.names = 1)
View(df)
# -----------


pdf(paste(output_dir,"barplot_TIC_befSum.pdf", sep = ""))
barplot(colSums(df, na.rm = TRUE), las=2)
dev.off()


pdf(paste(output_dir,"barplot_TIC_befSum_log2scale.pdf", sep = ""))
barplot(colSums(log(df,2), na.rm = TRUE), las=2)
dev.off()


pdf(paste(output_dir,"corr_plot_befSum_log2.pdf", sep = ""))
corr <- cor(as.matrix(log(df,2)), method = "pearson", use = "complete.obs")
corrplot(corr, order = "hclust", col = rev(brewer.pal(8,"RdYlBu")), type = "upper", cl.lim = c(0, 1),
         tl.col = "black", tl.srt = 90) # Text label color and rotation
dev.off()


# 
name_l <- sapply(colnames(df), function(x)unlist(strsplit(x,".mzML"))[1])
name_group_l <- c(name_l[1:6], sapply(name_l[7:length(name_l)], function(x)substring(x,1,nchar(x)-1))) # it is somekind of dictionary

group_iterator <- unique(name_group_l)

new_df <- as.data.frame(row.names(df))
for (item in group_iterator) {
  temp_df <- as.data.frame(df[,name_group_l %in% item]) #R automatically convert single column into non-df, so you need to convert it to df.
  if (nrow(temp_df) > 1) {
    new_col <- as.data.frame(rowMeans(temp_df, na.rm = TRUE))
  }
  else {
    new_col <- temp_df
  }
  new_df <- cbind(new_df,new_col)
}

colnames(new_df) <- c("FT_ID",group_iterator)
write.csv(new_df, paste(output_dir,"featureValues_summarized.csv"), row.names = FALSE)

#log2 scale
log2_new_df <- log(new_df[,2:ncol(new_df)],2)

write.csv(log2_new_df, paste(output_dir,"featureValues_summarized_log2scale.csv"))

#post-average bar plots of TIC
pdf(paste(output_dir,"barplot_TIC_postSum_log2.pdf", sep = ""))
barplot(colSums(log2_new_df, na.rm = TRUE), las=2)
dev.off()


pdf(paste(output_dir,"corr_plot_postSum_log2.pdf", sep = ""))
corr <- cor(as.matrix(log2_new_df), method = "pearson", use = "complete.obs") # new_df
corrplot(corr, order = "hclust", col = rev(brewer.pal(8,"RdYlBu")), type = "upper", cl.lim = c(0, 1),
         tl.col = "black", tl.srt = 90) # Text label color and rotation
dev.off()



# a principal componenet analysis
pc <- prcomp(t(na.omit(as.matrix(log2_new_df))), center = TRUE)

## Plot the PCA
group_colors <- c("#0000ff","#999999","#ff0000")
names(group_colors) <- c("G2_R5pos","G1_Naive","G3_R5neg")

cols <- group_colors[sapply(group_iterator,function(x)substring(x,1,nchar(x)-2))]
sample_name <- group_iterator

pdf(paste(output_dir,"PCA_with_labels_postSum_log2.pdf", sep = ""))
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2],  col = "darkgrey", labels = sample_name , #labels = "",
     pos = 1, cex = 1)
dev.off()

