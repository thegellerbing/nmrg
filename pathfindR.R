#Batch pathfindR Processing
library(pacman)
p_load(pathfindR, tidyverse, magrittr, readxl, writexl, org.Hs.eg.db)
path <- #Folder Path
files <- list.files(pattern = "_.xlsx") #Differential Expression File with HGNC Symbol

for (i in 1:length(files)){
  filename <- files[i]
  n <- gsub("_.xlsx", "", filename)
  cat("Now performing pathfindR for", n, "\n")
  df <- read_excel(filename)
  df %>% filter(padj != "NA") %>% as.data.frame() -> fildf
  fildf %<>% filter(hgnc_symbol != "NA")
  fildf %<>% dplyr::select(hgnc_symbol, log2FoldChange, padj)
  
  cat("Now running GO analysis \n")
  outgo <- run_pathfindR(fildf, gene_sets = "GO-BP", iterations = 100)
  cat("Now running KEGG analysis \n")
  outkegg <- run_pathfindR(fildf, gene_sets = "KEGG", iterations = 100)
  cat("Now running Reactome analysis \n")
  outrtm <- run_pathfindR(fildf, gene_sets = "Reactome", iterations = 100)
  
  assign(paste("go", n, sep = ""), outgo)
  assign(paste("kegg", n, sep = ""), outkegg)
  assign(paste("rtm", n, sep = ""), outrtm)
  
  name1 <- paste(path, n, "_go_siggenes.xlsx", sep = "")
  name2 <- paste(path, n, "_kegg_siggenes.xlsx", sep = "")
  name3 <- paste(path, n, "_reactome_siggenes.xlsx", sep = "")
  
  write_xlsx(outgo, name1)
  write_xlsx(outkegg, name2)
  write_xlsx(outrtm, name3)
}
rm(files, i, filename, df, fildf, outgo, outkegg, outrtm, name1, name2, name3, n, path)

#Batch pathfindR Diagram Processing
library(pacman)
p_load(pathfindR, tidyverse, magrittr, readxl, writexl)
path <- #Folder Path
files <- list.files(path = "", pattern = "_xlsx")
dataframe <- ls(pattern = "go")

for (i in 1:length(files)){
  filename <- files[i]
  n <- gsub("_shrinked_fullres.xlsx", "", filename)
  cat("Now processing", n, "\n")
  readfn <- paste("C:/Users/PC/OneDrive - Monash University/Desktop/PhD Folder/PhD Project/Data/RNA-Seq/Self-Aligned/DESeq Results/Shrinked LFC/apeglm/Filtered txt Files/", filename, sep = "")
  df <- read_excel(readfn)
  df %>% filter(padj != "NA") %>% as.data.frame() -> fildf
  fildf %<>% filter(hgnc_symbol != "NA")
  fildf %<>% dplyr::select(hgnc_symbol, log2FoldChange, padj)
  x <- dataframe[i]
  output <- get(x)
  
  visualize_terms(result_df = output, input_processed = fildf, hsa_KEGG = FALSE)
  file.rename("term_visualizations", n)
}
rm(files, dataframe, i, filename, readfn, df, fildf, x, output, n)