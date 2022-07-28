## Density Distribution

#load libraries 
shhh <- suppressPackageStartupMessages # It's a library, so shhh!
shhh(library(FactoMineR))
shhh(library(factoextra))
shhh(library(reshape2))
shhh(library(tidyverse))
shhh(library(RColorBrewer))
shhh(library(gridExtra))
shhh(library(pheatmap))
shhh(library(ggbiplot))
library(repr)        # jupyter notebook R ; image manipulation 
#' shhh(library(devtools))
#' install_github("vqv/ggbiplot")

options(warn=-1)     # don't show warnings in Output 
options(digits = 14) # Or any higher number; fixing decimal places

# Custom density distribution plot
#' this plot is aimed to plot density distribution of all features.
#' It provides opportunity to check for samples that are distinctively different (e.g., samples that have generally lower signals)

custom_density_plot <- function(data, y_lim_param) {
    #' data : - Data can be log2(preferably) /linear; if have NAs will be imputed to calculate density
    #' y_lim_param; to adjust height of y axis(density) in the plot

    #' Imputing for 0s
    data[is.na(data)] <- 0
    options(digits = 4)
    long <- melt(data[, 2:ncol(data)]) # remove the feature identifier information. 

    if (y_lim_param > 0) {

        ggplot(aes(x = value, colour = variable), # this color variable right here is imp to plot each sample individually
            data = long
        ) +
            geom_density() +
            theme_bw() +
            theme(legend.position = "") +
            scale_colour_grey() + # remove it for more colorful plot
            #' main title settings
            theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5)) +
            labs(title = "Density Distribution of Data") +
            #' axis settings
            theme(axis.text = element_text(size = 18)) +

            #' axis labels settings
            theme(axis.title = element_text(size = 20)) +

            #' axis labels
            xlab("Log2 Intensity") +
            ylab("Density") +
            coord_cartesian(xlim = c(7, 30), ylim = c(0, y_lim_param)) # adjust x-asis limits
    } else {

       ggplot(aes(x = value, colour = variable), # this color variable right here is imp to plot each sample individually
            data = long
        ) +
            geom_density() +
            theme_bw() +
            theme(legend.position = "") +
            scale_colour_grey() + # remove it for more colorful plot
            #' main title settings
            theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5)) +
            labs(title = "Density Distribution of Data") +
            #' axis settings
            theme(axis.text = element_text(size = 18)) +

            #' axis labels settings
            theme(axis.title = element_text(size = 20)) +

            #' axis labels
            xlab("Log2 Intensity") +
            ylab("Density") +
            coord_cartesian(xlim = c(7, 30)) # adjust x-asis limits
    }
}


# Plot TIC

plotTIC <- function(data, metadata, metadataSampleIdentifier, metaDataCol4color) {

    # data; log2/linear scale data with nas
    # metadata is optional
    # data[is.na(data)]<-0

    if (nrow(metadata) > 0) { # confusing, why nrow > 0, then it is plot without factors?
        # plot without factors


        df.tic <- data.frame(
            columnsum = colSums(data[, 2:ncol(data)], na.rm = TRUE),
            Sample.ID = names(data[, 2:ncol(data)])
        )

        df.tic <- merge(df.tic, metadata, by.x = "Sample.ID", by.y = metadataSampleIdentifier, all.x = TRUE)
        df.tic[, metaDataCol4color] <- as.factor(df.tic[, metaDataCol4color])

        ggplot(data = df.tic, aes(x = Sample.ID, y = columnsum, color = metaDataCol4color)) + #
            geom_bar(stat = "identity") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) + # rotate the axis
            theme(text = element_text(size = 10)) + # font size
            ggtitle("Mean Intensity") +
            xlab("FileName ") +
            ylab("TIC log2 Intensity") +
            theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 22),
                axis.title.y = element_text(size = 14)
            )
    } else {

        # plot with factors

        df.tic <- data.frame(
            columnsum = colSums(data[, 2:ncol(data)], na.rm = TRUE),
            Sample.ID = names(data[, 2:ncol(data)])
        ) #

        df.tic <- merge(df.tic, metadata, by.x = "Sample.ID", by.y = metadataSampleIdentifier, all.x = TRUE)
        ggplot(data = df.tic, aes(x = Sample.ID, y = columnsum)) + #
            geom_bar(stat = "identity") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) + # rotate the axis
            theme(text = element_text(size = 10)) + # font size
            ggtitle("Mean Intensity") +
            xlab("FileName ") +
            ylab("TIC log2 Intensity") +
            theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 14)
            )
    }
} 


# Get TIC outlier

#' `generate_TIC_table` creates a table that has `columnsum` & `Sample.ID` as columns
#'@Example 
#'	               columnsum	Sample.ID
#' SAM000203_tp1	20301	SAM000203_tp1
#' SAM000440_tp1	20334	SAM000440_tp1
#' SAM000427_tp1	20349	SAM000427_tp1


generate_TIC_table <- function(data) {
    # sdcutoff = standard deviation cutoff

    df.TIC <- data.frame(
        columnsum = colSums(data[, 2:ncol(data)], na.rm = TRUE),
        Sample.ID = names(data[, 2:ncol(data)])
    )


    return(df.TIC)
}

#' find outliers based how much of standard deviation you can tolerate from the mean of the overall dataset.
get_TIC_outliers <- function(data, Standard_deviation_count) {
    df.TIC <- generate_TIC_table(data)
    Outliers <- df.TIC %>% dplyr::filter(columnsum <  mean(df.TIC$columnsum) -  sd(df.TIC$columnsum)*Standard_deviation_count)
    View(Outliers)
}



# Plot Missing Values across sample
#' Plot number of missing data for a given feature.

plot_NA_cdf <- function(data, present_percentage = FALSE) {
    #' data is feature table with first column as identifier and rest every column is a sample
    #' logged 2/linear scale  data with NAs ; no imputation happening
    #' Calculate missing values count for each features
    #' Then create cumulative distribution of missing values.
    #' Options to present percentage of NA instead of NA count
    
    na_count <- apply(data[, 2:ncol(data)], 1, function(x) sum(is.na(x))) %>% data.frame()
    colnames(na_count) = c("na_count")
    if (present_percentage == FALSE) {
        ggplot(data = na_count, aes(na_count)) + stat_ecdf(geom = "point") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(text = element_text(size = 10)) +
        theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5)) +
        labs(title = "Missing Values Distribution") +
        # axis settings
        theme(axis.text = element_text(size = 16)) +

        # axis labels settings
        theme(axis.title = element_text(size = 20)) +
        xlab("Number of missing values") +
        ylab("Percentage of features") 
    } else {
        na_count$percentageNA = na_count$na_count/nrow(data)
        
        ggplot(data = na_count, aes(percentageNA)) + stat_ecdf(geom = "point") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(text = element_text(size = 10)) +
        theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5)) +
        labs(title = "Missing Values Distribution") +
        # axis settings
        theme(axis.text = element_text(size = 16)) +

        # axis labels settings
        theme(axis.title = element_text(size = 20)) +
        xlab(paste0("Percentage of NA per feature (Total: ", nrow(data),")")) +
        ylab("Percentage of features") 
         
    }
    
}



# Heatmap with annotation
#' plot heatmap such as correlation matrix

plot_heatmap <- function(data, metadata, metaDataIdentifier, metaDataCol4color) {

    #' data; log2/linear scale data with NAs
    #' metadata is optional
    #' data[is.na(data)]<- 0 

    if (nrow(metadata) > 0) {
        col_annot <- data.frame(FileName = names(data[, 2:ncol(data)])) # no mz and rt columns
        col_annot <- merge(col_annot, metadata, by.x = "FileName", by.y = metaDataIdentifier, x.all = TRUE)
        row.names(col_annot) <- col_annot$FileName
        col_annot <- col_annot %>% dplyr::select(metaDataCol4color)

        cormat <- cor(data[, 2:ncol(data)], use = "pairwise.complete.obs")
        # Sets the minimum (0), the maximum (1), and the increasing steps (+0.1) for the color scale
        # Note: if some of your data are outside of this range, they will appear white on the heatmap
        breaksList <- seq(0, 1, by = 0.1)


        pheatmap(cormat,
                cluster_row = FALSE, 
                cluster_col = FALSE, 
                fontsize = 15, 
                fontsize_row = 2,
                fontsize_col = 6, 
                scale = "none", 
                show_rownames = F, 
                show_colnames = T, 
                main = "Sample-Sample-Correlation",
                # cellwidth = 2, cellheight = 2
                treeheight_row = 40, 
                treeheight_col = 40, 
                annotation_col = col_annot,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList
        )
    } else {
        cormat <- cor(data[, 2:ncol(data)], use = "pairwise.complete.obs")
        breaksList <- seq(0, 1, by = 0.1)
        pheatmap(cormat,
                 cluster_row = FALSE, 
                 cluster_col = FALSE,
                 fontsize = 15, 
                 fontsize_row = 2, 
                 fontsize_col = 6, 
                 scale = "none",
                 show_rownames = F, 
                 show_colnames = T, 
                 main = "Sample-Sample-Correlation",
               # cellwidth=2, cellheight=2
                 treeheight_row = 40, treeheight_col = 40,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList
        )
    }
}


# Get outliers heatmap
#' 
#'
#'

get_outliers_heatmap <- function(data) {
    #'' - data; log2/linear scale data with nas
    col_annot <- data.frame(FileName = names(data[, 2:ncol(data)])) # no mz and rt columns

    # assuming first two columns are rt and mz
    cormat <- cor(data[, 2:ncol(data)], use = "pairwise.complete.obs")
    df.outlier <- data.frame(samples = colnames(cormat), columnmeancor = colMeans(cormat))

    # mean.df.outlier <-mean(df.outlier$columnmeancor)

    # sd.df.outlier<-sd(df.outlier$columnmeancor)

    # cutoff<-mean.df.outlier- sd.df.outlier
    cutoff <- 0.7
    df.outlier <- df.outlier %>% mutate(SampleType = ifelse(df.outlier$columnmeancor > cutoff, "inlier", "outlier"))
    df.outlier

    return(df.outlier)
}


getoutliersplot <- function(data) {
    ggplot(df.outliers, aes(y = columnmeancor, x = samples, colour = SampleType)) +
        geom_point(size = 4) +
        theme_bw() +
        xlab("") +
        ylab("Mean sample correlation") +
        theme(text = element_text(size = 8)) +
        theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5)) +
        labs(title = "Outliers Distribution") +
        # axis settings
        theme(axis.text.x = element_text(size = 8)) +
        # axis labels settings
        theme(
            axis.title = element_text(size = 12),
            axis.text.x = element_text(size = 10, angle = 90, face = "plain"),
            axis.text.y = element_text(size = 14, angle = 90, face = "plain")
        )
}


# Plot PCA

options(repr.plot.width = 18, repr.plot.height = 16, res = 200)
plotPCA <- function(data,
                    df.metadata = "nodata",
                    identifier,
                    columns4color,
                    contain_NA = FALSE) {

    ## pass  logged2 data with Nas ;
    # Missing values are imputed by the mean of the variable: you should use the imputePCA function of the missMDA package
    # default, contain NA = FALSE; If it is TRUE, then the current function replace NA with 0.

    if (df.metadata == "nodata") {
        print("No metadata; Try GAin with Metadata")
    } else {

        # Missing values are imputed by the mean of the variable: you should use the imputePCA function of the missMDA package
        # dim(df.featuretable)
        df.featuretable.pca <- t(data[, 2:ncol(data)]) %>% data.frame()
        colnames(df.featuretable.pca) <- data[, 1]
        # head(df.featuretable)

        df.featuretable.pca[[identifier]] <- row.names(df.featuretable.pca)

        colnames(df.metadata)

        df.metadata <- df.metadata %>% dplyr::select(identifier, columns4color) # `
        df.featuretable.pca <- left_join(df.featuretable.pca, df.metadata)

        ncol(df.featuretable.pca)
        tail(colnames(df.featuretable.pca))
        # head(colnames(df.featuretable.pca))

        if (contain_NA == TRUE) {
            df.featuretable.pca[is.na(df.featuretable.pca)] <- 0
        }
        # return(df.featuretable.pca)

        res.pca <- PCA(df.featuretable.pca[, 1:(ncol(df.featuretable.pca) - 2)], graph = FALSE)
        fviz_pca_ind(res.pca,
            geom.ind = "point", # show points only (nbut not "text")
            col.ind = df.featuretable.pca[[columns4color]], # color by groups
            palette = "jco", repel = TRUE,
            addEllipses = TRUE, # Concentration ellipses
            legend.title = "Groups"
        )
    }
}
