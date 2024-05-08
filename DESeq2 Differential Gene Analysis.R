#Installing Packages
install.packages('pacman')
install.packages('BiocManager')
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install("EnhancedVolcano")

#Loading Pacakages
library(pacman)
p_load(dplyr,tidyverse, DESeq2, data.table,survival, survminer, readxl, magrittr, biomaRt, writexl, apeglm, ggplot2, vsn, EnhancedVolcano, ggbeeswarm, RColorBrewer, pheatmap, ashr, BiocParallel)
register(SnowParam(4))

#Sourcing Functions
source("C:/Users/PC/OneDrive - Monash University/Desktop/PhD Folder/Projects/R/Codes/Functions/deseq_func.R") 

#Function Description
#tidyexpr helps convert log2 transformed count data into raw counts
#tidymd produces a metadata file required for the differential expression
#rundeseq runs the differential expression analysis
#getresults processes the differential expression analysis and allows you to save them into excel files
#shrinklfc shrinks the log fold change of the differential expression analysis (an extra layer of normalisation)
#getlfcresults processes the shrinked differential expression analysis and allows you to save them into excel files
#maplot allows you to plot individual MA plots
#volcanoplot allows you to plot individual volcano plots
#countplot allows you to plot individual count plot for specific genes
#datatransform further transforms the differential expression analysis for further downstream analysis such as PCA plots and sample to sample distance heatmaps
#cheatmap allows you to plot expression heatmap of your gene of interest

## expr -> normalised count data
## md -> produced metadata from the tidyexpr function


#If Data Transformation is required
expr <- tidyexpr(expr)

#To establish the clinical data required for the differential expression analysis
md <- tidymd(expr)

#Running Differential Expression Analysis
dds <- rundeseq(expr, md)
ensembl = useMart( "ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org", dataset="hsapiens_gene_ensembl" )
attrlist = listAttributes(ensembl) #Run this and view the list for the appropriate tag if the genes in your provided expression data aren't in official HUGO symbols or Ensembl gene IDs
results <- getresults(dds)
lfc <- shrinklfc(dds)
lfcres <- getlfcresults(lfc)

#MA Plot (Please provide the results from 'getresults'/'shrinklfc')
maplot(results/lfc)

#Volcano Plot
volcanoplot(results)

#Count Plot
countplot(dds)

#Data Transformation and Visualisation
tres <- datatransform(dds)

#Expression Heatmap of Specific Genes
cheatmap(tres)

genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = c("ensembl_gene_id"),
                 values = results$ensembl,
                 mart = ensembl)