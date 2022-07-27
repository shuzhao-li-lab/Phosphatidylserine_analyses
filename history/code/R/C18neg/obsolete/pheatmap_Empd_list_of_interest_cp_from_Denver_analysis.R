epd_l <- read.table("../../output/RPneg/RPneg_padj_Denver/tables/epd_of_interest.txt", header = 1)

lofCpd_df <- read.table("../../output/RPneg/RPneg_padj_Denver/tables/ListOfEmpiricalCompounds_mod.txt", 
                        header = 1,fill = TRUE)
featab <- read.csv("../../output/RPneg/featureValues.csv")


m_df <- merge(epd_l,lofCpd_df, by.x = "epd_list", by.y = "EID", how = "inner")
m_df$FT_index <- gsub("row", "",m_df$massfeature_rows)
FTindex = unique(gsub("row", "",m_df$massfeature_rows))

m_df <- merge(featab, m_df, by.x = 0, by.y = "FT_index")

m_df$cbn_name <- paste0(m_df$str_row_ion,m_df$compound_names,sep = "|")
row.names(m_df) <- m_df$cbn_name
colnames(m_df)
HT <- m_df[,3:42]

log2HT <- log(HT,2)

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
annotation_col = data.frame(label = factor(c(rep("HEU",20),rep("HUU",20)))) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(z_transformed)  #e.g., PUF1_1, PUF_2, Dpmt1

#annotation_row = data.frame(GeneClass = factor(HT$pathways)) # metabolism A...
#rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
label_color = c("#ff0000","#0000ff")
names(label_color) = c("HEU","HUU") #!
label_color
# HEU       HUU 
# "#ff0000" "#0000ff" 

range(z_transformed) #check the range of z_transformed
breaksList = seq(-2.4, 2.4, 0.1)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()


mat_breaks <- seq(min(z_transformed), max(z_transformed), length.out = 1)

out_pdf_file = paste("mcg_RPneg_padj05",".pdf",sep = "")  #change the number


pdf(file=out_pdf_file, width =15,height=10, paper = "special",onefile=FALSE)
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = NULL, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()