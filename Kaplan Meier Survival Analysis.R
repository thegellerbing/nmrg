#Installing Packages
install.packages('pacman')
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('survival')
p_load(tidyverse, data.table, survminer, magrittr, readxl, writexl, maxstat)

#Loading And Sourcing Packages 
library(pacman)
p_load(tidyverse, data.table, survival, survminer, magrittr, readxl, writexl, maxstat, dplyr)


source("C:/Users/PC/OneDrive/Desktop/PhD Folder/Projects/R/Codes/Functions/KM_func.R") 

#Function Description
#genes2hgnc converts your expression data into official HUGO symbol if desired
#tidydata prepares a dataframe with all the relevant data required for the survival analysis
#umiKM runs the univariate survival analysis, it returns plots and a table/list and saves the results into individual excel files
#multicox runs the multivariate survival analysis

exprdata #Normalised Expression Data
md #Survival Data

#genes2hgnc function converts ensembl IDs to official gene symbols. 
ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org", dataset="hsapiens_gene_ensembl" )
exprdata <- genes2hgnc(exprdata) 
rm(ensembl)

#tidydata function prepares the dataframe that is required for the survival analyses
#Requires a normalised expression dataframe, a clinical information dataframe and a dataframe containing a list of genes
kmdf <- tidydata()

#Univariate KM Analysis####
unidf <- uniKM(kmdf)

#Multivariate Cox Analysis ####
cox <- multicox(kmdf)
