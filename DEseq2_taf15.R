## DESeq2 on FINAL taf15 data. 
## Count files were generated using htseq-count on genepool.

## The project now contains data from UC and MO from stage 10 and 15. 
## each set done in triplicate except UC stage 15 which has four replicates. 

## Load libraries.
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")

library( DESeq2 ) #1.6.3, up-to-date
library( gplots )
library( ggplot2 )
library( RColorBrewer )
sessionInfo()

## Create palettes
pccol <- colorRampPalette( brewer.pal( 10, "Spectral" ) ) ( 6 ) # for PCA plots. I should make this less candy.
distcol <- colorRampPalette( rev( brewer.pal( 9, "GnBu" ) ) ) (20) # For distance plots
hmcolMono20 <- colorRampPalette ( brewer.pal( 9, "GnBu") ) ( 20 ) # for one color heatmaps
hmcolBR <- colorRampPalette ( rev( brewer.pal( 9, "RdBu" ) ) ) ( 20 ) # For two color heatmaps

## Set local working directory
setwd("./")

## Create design matrix
inDir <- file.path( paste( getwd(), "HTSeqCount_TAF15/", sep="/" ) )
countFiles <- list.files( path= inDir, pattern= "(st10_UC*)|(st10_MO*)|(st15_UC*)|(st15_MO*)" )
sampleIDs <- c( "UC1_ST10", "UC2_ST10", "UC3_ST10", 
                "MO1_ST10", "MO2_ST10", "MO3_ST10",
                "UC1_ST15", "UC2_ST15", "UC3_ST15", "UC4_ST15",
                "MO1_ST15", "MO2_ST15", "MO3_ST15" )
conditions <- c( rep( "CTRL_ST10", 3 ), 
                 rep( "TAF15MO_ST10", 3 ), 
                 rep( "CTRL_ST15", 4 ), #Remember the bonus sample
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
                                     #design= ~condition+stage )
                                     design= ~condition )

### Sanity plots
### First the transformations
ddsTaf <- DESeq( ddsTaf )        
rldTaf <- rlog( ddsTaf )
vsdTaf <- varianceStabilizingTransformation( ddsTaf )

### PCA 
### BROKEN
plotPCA( rldTaf, intgroup= c("condition", "stage") ) + theme_bw() + geom_text( aes( label=names ), hjust= 0.5, vjust= -0.75 ) + geom_point( size= 3 )

### Distance plots
distTaf <- dist( t( assay( rldTaf ) ) ) 
distTaf.mat <- as.matrix( distTaf )
heatmap.2 (distTaf.mat, trace="none", col= distcol, Colv=T, Rowv=F, main= "Distance Matrix of TAF15 Experiment" )


###################
## Getting the results
## Prepared separate sets for all DE genes and DE genes +2FC

###################
## For stage 10
res10 <- results( ddsTaf, contrast= c( "condition", "CTRL_ST10", "TAF15MO_ST10" ) )
table(  res10$padj < 0.1 ) # 169 genes, surprisingly few
table( res10$padj < 0.1 & res10$log2FoldChange > abs( 1 ) ) # 34 genes
sig10 <- subset( res10, padj < 0.1 )
sig10.2fc <- subset( res10, padj < 0.1, log2FoldChange > abs( 1 ) )

###################
## For stage 15
res15 <- results ( ddsTaf, contrast= c( "condition", "CTRL_ST15", "TAF15MO_ST15" ) )
table( res15$padj < 0.1 ) # 1853 surprisingly many
table( res15$padj < 0.1 & res15$log2FoldChange > abs( 1 ) ) #544
sig15 <- subset( res15, padj < 0.1 )
sig15.2fc <- subset( res15, padj < 0.1 & log2FoldChange > abs( 1 ) )

### Intersect DE genes from stage 10 and 15
## Find common genes
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

### Write files for Caitlin to look at.
write.table( x = as.data.frame( sig10 ), file= "DEGenes4Caitlin/sigTaf15_stage10.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig10.2fc ), file= "DEGenes4Caitlin/sigTaf15_stage10_2FC.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig15 ), file= "DEGenes4Caitlin/sigTaf15_stage15.txt", sep= "\t", quote= F )
write.table( x= as.data.frame( sig15.2fc ), file= "DEGenes4Caitlin/sigTaf15_stage15_2FC.txt", sep= "\t", quote= F )
write.table( x = intersectDEG, file= "DEGenes4Caitlin/intersect_DEGenes.txt", sep= "\t", quote= F )

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
