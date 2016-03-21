## DESeq2 on FINAL taf15 data. 
## Count files were generated using htseq-count on genepool.

## The project now contains data from UC and MO from stage 10 and 15. 
## each set done in triplicate except UC stage 15 which has four replicates. 


#####################################
## LOAD REQUIRED LIBRARIES
#####################################

### Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
library( DESeq2 ) # 1.10.1

### Universal packages
library( gplots ) # 2.17.0
library( ggplot2 ) # 2.1.0
library( RColorBrewer ) # 1.1-2

sessionInfo()



### TODO: Move palettes to relevant sections

### Create heatmap color palette
hmcol <- colorRampPalette( rev( brewer.pal( 9, "RdBu" ) ) ) ( 255 )
### Create distance matrix color palette
distcol <- colorRampPalette( rev( brewer.pal(9, "GnBu") ) ) ( 20 )

## Set local working directory
# setwd("./")

### Create design matrix
# inDir <- file.path( paste( getwd(), "HTSeqCount_TAF15/", sep="/" ) )
inDir <- file.path( "HTSeqCount_TAF15/" )
countFiles <- list.files( path= inDir, pattern= "(st10_CTR*)|(st10_MO*)|(st15_CTR*)|(st15_MO*)" )
# sampleIDs <- c( "UC1_ST10", "UC2_ST10", "UC3_ST10", 
#                 "MO1_ST10", "MO2_ST10", "MO3_ST10",
#                 "UC1_ST15", "UC2_ST15", "UC3_ST15", "UC4_ST15",
#                 "MO1_ST15", "MO2_ST15", "MO3_ST15" )
sampleIDs <- c( "CTR1_ST10", "CTR2_ST10", "CTR3_ST10", 
                "MO1_ST10", "MO2_ST10", "MO3_ST10",
                "CTR1_ST15", "CTR2_ST15", "CTR3_ST15", "CTR4_ST15",
                "MO1_ST15", "MO2_ST15", "MO3_ST15" )

conditions <- c( rep( "CTRL_ST10", 3 ), 
                 rep( "TAF15MO_ST10", 3 ), 
                 rep( "CTRL_ST15", 4 ), # Remember the extra sample
                 rep( "TAF15MO_ST15", 3 ) )
stages <- c (rep( "stage10", 6 ), 
           rep( "stage15", 7 ) )
tafTable <- data.frame( sampleName= sampleIDs, 
                          fileName= countFiles, 
                          condition= conditions, 
                          stage= stages )

tafTable # Looks good
rm( sampleIDs, conditions, stages ) # Clean up.


ddsTaf <- DESeqDataSetFromHTSeqCount( sampleTable = tafTable, 
                                     directory = inDir, 
                                     design= ~condition )


### First the transformations
ddsTaf <- DESeq( ddsTaf )        
rldTaf <- rlog( ddsTaf )
vsdTaf <- varianceStabilizingTransformation( ddsTaf )

############################
###### Sanity plots ########
############################

### Create directory for plots
dir.create("plots/")

### PCA Plot
### Colors from http://colorbrewer2.org (4-class RdYlBu)
### Note: the colors are assigned to categories alphabetically.
pccol <- c( '#d7191c', '#2c7bb6', '#fdae61', '#abd9e9' )
pc <- plotPCA( rldTaf, returnData = TRUE )

pca <- ggplot(data = pc, aes( pc[,1], pc[,2] ) ) + theme_bw()
pca <- pca + geom_point( size = 4, alpha = 0.7, aes( colour = group ) )   
pca <- pca + scale_colour_manual(values = pccol, guide_legend("") ) 
pca <- pca + ggtitle( "PCA plot\nrld transformed values")
pca <- pca + labs(x = "PC1", y = "PC2" )
pca <- pca + theme(legend.position = c(0.15, 0.9))
# pca
ggsave("plots/PCA.pdf", width = 8, height = 8)


### Distance plot
### Create distance matrix color palette
### TODO: Make this pretty. Consult with somite code.
distcol <- colorRampPalette( rev( brewer.pal(9, "GnBu") ) ) ( 10 )

distTaf <- dist( t( assay( rldTaf ) ) ) 
distTaf.mat <- as.matrix( distTaf )
heatmap.2( distTaf.mat, trace="none", 
           col= distcol, 
           Colv=FALSE, Rowv=FALSE, 
           main= "Distance Matrix of TAF15 Experiment" )


################################################
### Get the results
### Prepared separate sets for all DE genes and DE genes +2FC
################################################

################################################
### Stage 10 results
### Contrast: condition, numerator, denominator.
################################################
res10 <- results( ddsTaf, contrast= c( "condition", "TAF15MO_ST10", "CTRL_ST10" ) )
table(  res10$padj < 0.1 ) # 296 genes
table( res10$padj < 0.1 & res10$log2FoldChange > abs( 1 ) ) # 82 genes
sig10 <- subset( res10, padj < 0.1 )
sig10.2fc <- subset( res10, padj < 0.1, log2FoldChange > abs( 1 ) )

################################################
### Stage 15 results
################################################
res15 <- results ( ddsTaf, contrast= c( "condition", "TAF15MO_ST15", "CTRL_ST15" ) )
table( res15$padj < 0.1 ) # 4352, surprisingly many
table( res15$padj < 0.1 & res15$log2FoldChange > abs( 1 ) ) # 1148 
sig15 <- subset( res15, padj < 0.1 )
sig15.2fc <- subset( res15, padj < 0.1 & log2FoldChange > abs( 1 ) )


### Identify shared DE genes from stage 10 and 15
commonGenes <- intersect( sig10@rownames, sig15@rownames )

## Make stage-specific dataframe for merging
df10 <- as.data.frame( sig10 )
df10 <- df10[ rownames( df10 ) %in% commonGenes, ]
df15 <- as.data.frame( sig15 )
df15 <- df15[ rownames( df15 ) %in% commonGenes, ]

## Merge data frames by row (gene) names
intersectDEG <- merge( df10, df15, by = "row.names", suffixes = c( ".st10", ".st15" ) )
rownames(intersectDEG) <- intersectDEG[, 1 ]
intersectDEG$Row.names <- NULL


###############################
### Write lists of DE genes
###############################

### Create directory for tables
dir.create("DE_genesList")

### Write the tables
write.table( x = as.data.frame( sig10 ), file= "DE_genesList/sigTaf15_stage10.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig10.2fc ), file= "DE_genesList/sigTaf15_stage10_2FC.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig15 ), file= "DE_genesList/sigTaf15_stage15.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig15.2fc ), file= "DE_genesList/sigTaf15_stage15_2FC.txt", sep= "\t", quote= F )
write.table( x = intersectDEG, file= "DE_genesList/intersect_DEGenes.txt", sep= "\t", quote= F )


############################################################
############################################################
############# TESTED CODE TO HERE 3/19/2016
############################################################
############################################################


#################
rm( df10, df15 ) # Clean up.
#################

#################
## Coming back to the analysis after two crazy househunting weeks
## Darwin 4 May, 2015
#################

#################
## Make heatmap of DE genes
#################

## First let's go back and test transformations
## Test the diffferent variance transformations
## Lowest variation with rld, so we go with that.
library("vsn") #vignette("vsn")
par(mfrow = c(1, 3))
notAllZero <- ( rowSums( counts( ddsTaf ) ) > 0 )
meanSdPlot( log2 (counts ( ddsTaf, normalized = T )[ notAllZero, ] +1 ), main="Log2 transformed")
meanSdPlot( assay( rldTaf[ notAllZero, ] ), main="Regularized log transformation" )
meanSdPlot( assay( vsdTaf[ notAllZero, ] ), main= "Variance Stabilizing transformation" )
par(mfrow= c(1, 1))

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
