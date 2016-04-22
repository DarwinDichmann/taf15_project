# DESeq2 on FINAL taf15 data. 
# Count files were generated using htseq-count on genepool.

# The project now contains data from UC and MO from stage 10 and 15. 
# each set done in triplicate except UC stage 15 which has four replicates. 



# Load libraries -----------------------

# Bioconductor packages
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
library(DESeq2) # 1.10.1

# Universal packages
library(gplots) # 2.17.0
library(ggplot2) # 2.1.0
library(RColorBrewer) # 1.1-2
library(vsn)  # 3.38.0
library(cowplot)  # 0.6.1
sessionInfo()


# TODO: Move palettes to relevant sections
# Create heatmap color palette
# hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu"))) (255)

### Function definitions --------------
writeDESeqResults <- function(dds, num, denom, 
                              condition = "condition", 
                              log2_cut = 1,  # 
                              padj_cut = 0.1,
                              return_data = TRUE,
                              filename = NULL) {
        # Convenience function that extracts results from a DESeqResult 
        # object and returns a dataframe or writes a file if a filename is
        # provided.
        #  
        # Args:
        #       dds:            A DESeqDataSet.
        #       num:            Group used as numerator in contrast.
        #       denom:          Group used denominator in contrast.
        #       condition:      Factor be used for contrast. Default: "condition".
        #       log2_cut:       Log2FC threshold for excluding genes. Default: 1.
        #       padj_cut:       Adjusted p-value threshold. Default: 0.1. 
        #       return_data:    If TRUE, returns results as dataframe.
        #       filename:       Results will be written to file if provided.
        #                       Can be combined with return_data = TRUE.
        condition <- as.character(condition)
        num <- as.character(num)
        denom <- as.character(denom)
        deg <- results(dds, c(condition, num, denom))
        deg <- as.data.frame(deg)
        deg <- subset(deg, padj < padj_cut & abs(log2FoldChange) > log2_cut)
        deg <- deg[order(- deg$log2FoldChange), ]
        deg <- cbind(Row.Names = rownames(deg), deg)
        colnames(deg) [1] <- "Gene"
        if (! is.null(filename)) {
                write.table(x = deg, 
                            file = filename, 
                            quote = FALSE, 
                            sep = "\t")
        }
        if (return_data == TRUE) {
                return(deg)
        }
}


### DESIGN MATRIX ---------------------------------
count_dir <- file.path("HTSeqCount_TAF15/")
count_files <- list.files( path = count_dir, 
                          pattern = "(st10_CTR*)|(st10_MO*)|(st15_CTR*)|(st15_MO*)" )

sampleIDs <- c("CTR1_ST10", "CTR2_ST10", "CTR3_ST10", 
               "MO1_ST10", "MO2_ST10", "MO3_ST10",
               "CTR1_ST15", "CTR2_ST15", "CTR3_ST15", "CTR4_ST15", 
               "MO1_ST15", "MO2_ST15", "MO3_ST15")

conditions <- c(rep("CTRL_ST10", 3), 
                rep("TAF15MO_ST10", 3), 
                rep("CTRL_ST15", 4), # Remember the extra sample
                rep("TAF15MO_ST15", 3))

stages <- c(rep("stage10", 6), 
            rep("stage15", 7))

taf_table <- data.frame(sampleName = sampleIDs, 
                       fileName   = count_files, 
                       condition  = conditions, 
                       stage      = stages)

# taf_table  # Looks good
rm(count_files, count_dir, sampleIDs, conditions, stages)  # Clean up.


# DESEQ AND TRANSFORMATIONS ---------------------
taf_dds <- DESeqDataSetFromHTSeqCount(sampleTable = taf_table, 
                                      directory   = inDir, 
                                      design      = ~ condition)
taf_dds <- DESeq(taf_dds)        
taf_rld <- rlog(taf_dds)
taf_vsd <- varianceStabilizingTransformation(taf_dds)


### Sanity plots ----------------------------------------------------------
dir.create("plots/")  # Manuscript figures go here.

# Determine which transformation minimizes variation.
not_all_zero <- rowSums(counts(taf_dds)) > 0
plot_log <- meanSdPlot(log2(counts(taf_dds, 
                                   normalized = TRUE)[not_all_zero, ] + 1),
                       plot = FALSE)
plot_rld <- meanSdPlot(assay(taf_rld[not_all_zero, ]), 
                       plot = FALSE)
plot_vsd <- meanSdPlot(assay(taf_vsd[not_all_zero, ]), 
                       plot = FALSE)
plot_sd <- plot_grid(plot_log$gg, plot_rld$gg, plot_vsd$gg, 
                     labels = c("A", "B", "C"), 
                     ncol = 3, align = 'h')  # similar variance.
plot_sd  # vsd and rlog variations are similar. I will use rlog. 
rm(plot_log, plot_rld, plot_vsd, plot_sd)  # Clean up.
# rm(not_all_zero)  # Might need this later.


### Principal Component Plot ----------------------------------------------
# Colors from http://colorbrewer2.org (4-class RdYlBu)
# Note: the colors are assigned to categories alphabetically.
# TODO: fix order on guide legend (consult with somite code)
pc_col <- c('#d7191c', '#2c7bb6', '#fdae61', '#abd9e9')
pc_data <- plotPCA(taf_vsd, returnData = TRUE)
pc_plot <- ggplot(data = pc_data, aes( pc_data[, 1], pc_data[, 2])) 
pc_plot <- pc_plot + theme_cowplot()  # white theme, no lines.
pc_plot <- pc_plot + geom_point(size = 4, alpha = 0.7, aes(color = group))
pc_plot <- pc_plot + scale_colour_manual(values = pc_col, guide_legend(""))
pc_plot <- pc_plot + ggtitle("PCA plot\nvsd transformed values")
pc_plot <- pc_plot + labs(x = "PC1", y = "PC2" )
pc_plot <- pc_plot + theme(legend.position = c(0.2, 0.8))
# pc_plot
ggsave("plots/PCA.pdf", width = 6, height = 4)
rm(pc_col, pc_data, pc_plot)


### Distance plot - tells same story as PC plot.---------------------------

# dist_col <- colorRampPalette(rev(brewer.pal(9, "GnBu"))) (20)
# dist_taf <- dist(t(assay(taf_rld))) 
# dist_taf <- as.matrix(dist_taf)
# heatmap.2(dist_taf, 
#           trace = "none", 
#           col   = dist_col, 
#           Colv  = FALSE, 
#           Rowv  = FALSE, 
#           main  = "Distance Matrix of TAF15 Experiment")
# rm(dist_taf, dist_col)  # clean up.


##### Write DEG results ---------------------------------------------------

dir.create("DEG_lists")

# Stage 10, only +2FC
deg_10_2fc <- writeDESeqResults(dds = taf_dds, 
                            num = "TAF15MO_ST10", 
                            denom = "CTRL_ST10", 
                            log2_cut = 1, 
                            return_data = TRUE, 
                            filename = "DEG_lists/DEG_2fc_taf15_st10.txt")
# Stage 10, all DEG
deg_10 <- writeDESeqResults(dds = taf_dds, 
                            num = "TAF15MO_ST10", 
                            denom = "CTRL_ST10", 
                            log2_cut = 0, 
                            return_data = TRUE, 
                            filename = "DEG_lists/DEG_taf15_st10.txt")

# Stage 15, only +2FC
deg_15_2fc <- writeDESeqResults(dds = taf_dds, 
                                num = "TAF15MO_ST15", 
                                denom = "CTRL_ST15", 
                                log2_cut = 1, 
                                return_data = TRUE, 
                                filename = "DEG_lists/DEG_2fc_taf15_st15.txt")

# Stage 15, all DEG
deg_15 <- writeDESeqResults(dds = taf_dds, 
                            num = "TAF15MO_ST15", 
                            denom = "CTRL_ST15", 
                            log2_cut = 0, 
                            return_data = TRUE, 
                            filename = "DEG_lists/DEG_taf15_st15.txt")


# TODO: Is this still necessary? =====
res_10 <- results(taf_dds, 
                  contrast = c("condition", "TAF15MO_ST10", "CTRL_ST10"))
table(res_10$padj < 0.1)  # 292 genes
table((res_10$padj < 0.1) & abs(res_10$log2FoldChange) >= 1)  # 135 genes
sig_10 <- subset(res_10, padj < 0.1)  
sig_10_2fc <- subset(res_10, padj < 0.1 & abs(log2FoldChange) >= 1)

# Create dfs for exporting DEG tables
deg_10 <- as.data.frame(sig_10)
deg_10 <- deg_10[order(- deg_10$log2FoldChange), ]
deg_10_2fc <- as.data.frame(sig_10_2fc)
deg_10_2fc <- deg_10_2fc[order(- deg_10_2fc$log2FoldChange), ]


res_15 <- results(taf_dds, 
                  contrast = c("condition", "TAF15MO_ST15", "CTRL_ST15"))
table(res_15$padj < 0.1) # 4354, surprisingly many
table(res_15$padj < 0.1 & abs(res_15$log2FoldChange) >= 1) # 2094 
sig_15 <- subset(res_15, padj < 0.1)
sig_15_2fc <- subset(res_15, padj < 0.1 & abs(log2FoldChange) >= 1)


##### TODO: Figure out if this section is necessary or not. ===============
### Identify shared DE genes from stage 10 and 15
deg_common <- intersect(sig_10@rownames, sig_15@rownames)

## Make stage-specific dataframe for merging

df10 <- df10[rownames(df10) %in% commonGenes, ]
df15 <- as.data.frame(sig15)
df15 <- df15[ rownames(df15) %in% commonGenes, ]

## Merge data frames by row (gene) names
intersectDEG <- merge(df10, df15, by = "row.names", 
                      suffixes = c(".st10", ".st15"))
rownames(intersectDEG) <- intersectDEG[, 1]
intersectDEG$Row.names <- NULL


###############################
### Write lists of DE genes
###############################

### Create directory for tables


### Write the tables
write.table(x = as.data.frame(sig10), 
            file = "DE_genesList/sigTaf15_stage10.txt", 
            sep = "\t", quote = FALSE)
write.table(x = as.data.frame(sig10_2fc), 
            file = "DE_genesList/sigTaf15_stage10_2FC.txt", 
            sep = "\t", quote = FALSE)
write.table(x = as.data.frame(sig15), 
            file = "DE_genesList/sigTaf15_stage15.txt", 
            sep = "\t", quote = FALSE)
write.table(x = as.data.frame(sig15_2fc), 
            file = "DE_genesList/sigTaf15_stage15_2FC.txt", 
            sep = "\t", quote = FALSE)
write.table(x = intersectDEG, 
            file = "DE_genesList/intersect_DEGenes.txt", 
            sep = "\t", quote = FALSE)


############################################################
############################################################
############# TESTED CODE TO HERE 3/19/2016
############################################################
############################################################


#################
rm(df10, df15) # Clean up.
#################

#################
## Coming back to the analysis after two crazy househunting weeks
## Darwin 4 May, 2015
#################

#################
## Make heatmap of DE genes
#################



##################
## Heatmaps of DE genes @ stage 10
##################

sig15.2fc_genes <- row.names(sig15.2fc)
sig15.2fc_genes

rldTaf.mat <- as.matrix( assay( rldTaf ) )
head(rldTaf.mat)
## Collapse the replicates
library(limma)
design <- cbind( CTR_ST10 = c( rep( 1, 3 ), rep( 0, 10) ), 
                 TAF15MO_ST10 = c( rep( 0, 3), rep( 1, 3 ), rep( 0, 7) ), 
                 CTR_ST15 = c( rep( 0, 6), rep( 1, 4), rep( 0, 3 ) ),
                 TAF15MO_ST15 = c( rep( 0, 10), rep( 1, 3 ) )
                 )
rldCollapse.mat <- as.matrix( lmFit( assay( rldTaf ), design = design ) )
rldCollapse.mat.st10 <- rldCollapse.mat[, 1:2 ]
rldCollapse.mat.st15 <- rldCollapse.mat[ , 3:4 ]
head (rldCollapse.mat.st10)

summary(rldCollapse.mat.st15)
plot( rldCollapse.mat.st15 )


heatmap.2(  rldCollapse.mat.st15[ sig15.2fc_genes, ], 
          trace="none", 
          scale = "none",
          dendrogram= "row", 
          Colv=F, 
          Rowv=T,
          mar= c(6, 10),
          main= "diff genes",
          col = hmcolMono20
          )




str(rldTaf)
