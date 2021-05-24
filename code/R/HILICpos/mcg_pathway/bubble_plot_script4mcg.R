#before any operations. Please remove all empty space for the rows and columns

# setwd("/Users/gongm/Documents/projects/Denver_HEU_Pilot/analysis_4_03212021wetRun/code/R")

output_dir <- "../../data/output/v0412_fix/HILICpos/mcg_HILICpos/R5posvsR5neg/1618781765.861877.mcg_HILICpos_padj_0.05_R5posvsR5neg/mcg_pathway_bubble_plot/"
dir.create(output_dir)
df_filt <- read.table("../../data/output/v0412_fix/HILICpos/mcg_HILICpos/R5posvsR5neg/1618781765.861877.mcg_HILICpos_padj_0.05_R5posvsR5neg/tables/mcg_pathwayanalysis_mcg_HILICpos_padj_0.05_R5posvsR5neg.tsv", 
                      quote = "\"",, sep = "\t",header = TRUE)
summary(df_filt)

df_filt$log10pvalue <- -log(df_filt$p.value,10)

if(TRUE) {
  df_filt <- df_filt[df_filt$p.value < 0.05,]
}
dim(df_filt)

if(FALSE) {
  df_filt$group <- factor(df_filt$group) 
} else {df_filt$group <- "R5posvsR5neg"}

if(FALSE) {
  df_filt$color <- factor(df_filt$color)
} else {df_filt$color <- "#000000"}



df_filt$pathway = ordered(df_filt$pathway, level = rev(df_filt$pathway))
range(df_filt$log10pvalue)
#[1] 1.326042 2.271997
range(df_filt$overlap_size)
#[1] -3.386694  2.951662


library(ggplot2)
 
    pdf(paste(output_dir, "ggplot_bubble_p0.05_module.pdf",sep = ""), width = 5.2, height =2.8)
    ggplot(df_filt, aes(x=group, y=pathway)) +
    geom_point(aes(size=overlap_size, colour = log10pvalue)) +
    scale_size_continuous(range = c(2, 6)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", breaks = seq(1,2),0.1) +
    # set transparency
    # https://ggplot2.tidyverse.org/reference/theme.html
    theme(
      panel.grid.major = element_line(colour = "grey50",linetype = "dashed", size = 0.2),
      panel.border = element_rect(fill = NA),
      #panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )
dev.off()
    