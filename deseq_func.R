tidyexpr <- function(expr){
  print("Transformation Ongoing...")
  suppressWarnings({
    gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    while(is.na(gene_col) == TRUE || gene_col > ncol(expr)){
      message("Error detected, please try again :)")
      gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    }
  })
  colnames(expr)[gene_col] <- "geneid"
  expr %<>% relocate(geneid)
  id <- expr[,1] %>% as.data.frame()
  colnames(id) <- "gene"
  fun <- function(x){2^x}
  expr <- apply(expr[, 2:ncol(expr)], 2, fun)
  pseudo <- as.numeric(readline(prompt = "Pseudocount value: "))
  pseudofun <- function(x){as.integer(x-pseudo)}
  expr <- apply(expr, 2, pseudofun)
  expr %<>% as.data.frame()
  expr <- cbind(id, expr)
  return(expr)
  cat("Done :>")
}

tidymd <- function(expr){
  message("This function creates a metadata table required for the differential expression analysis")
  Sys.sleep(0.5)
  message("If a metadata is already present with the experimental conditions, feel free to skip this function \n")
  Sys.sleep(0.5)
  cat("This function runs on the basis of common identifiers that R can search for for each experimental condition\n")
  Sys.sleep(0.5)
  cat("\nFor example the dataframe below contains 2 experimental conditions and one control: \n")
  
  example_df <- data.frame(id = c("TSPAN6", "TNMD", "DPM1", "SCYL3"), 
                        c1 = c(479, 326, 792, 502), 
                        c2 = c(325, 12, 1148, 651), 
                        c3 = c(286, 10, 748, 233), 
                        c4 = c(422, 4, 614, 336),
                        c5 = c(327, 12, 536, 383)) 
  colnames(example_df)[2:ncol(example_df)] <- c("GTEX-WZTO-2926-SM", "GTEX-PVOW-2526-SM", "TCGA-13NYS-3126-SM", "TCGA-11NUK-2926-SM", "lgg-11DXW-1126-SM")
  print(example_df)
  Sys.sleep(0.5)
  message("The example dataframe has one experimental control and two experimental conditions, GBM and LGG")
  Sys.sleep(0.5)
  cat("\nThe metadata table produced from this function for the example would be as below:\n")
  example_df2 <- data.frame(id = c("GTEX-WZTO-2926-SM", "GTEX-PVOW-2526-SM", "TCGA-13NYS-3126-SM", "TCGA-11NUK-2926-SM", "lgg-11DXW-1126-SM"),
                            condition = c("Control", "Control", "GBM", "GBM", "LGG"))
  print(example_df2)
  Sys.sleep(0.5)
  cat("\nTo produce the file, you are required to provide unique identifier for each of your experimental condition\n")
  cat("Referring to the example, the unique identifier for the control samples would be 'GTEX'\n")
  cat("Similarly, the unique identifier for the GBM samples would be 'TCGA'\n")
  cat("Lastly, the unique identifier for the LGG sample would be 'lgg'\n")
  message("Please make sure the column names of your samples in the count data has unique identifiers between different experimental conditions")
  
  proc <- readline("Hopefully all is clear, if so, do you want to proceed? (yes/no): ")
  
  if (tolower(proc) == "yes" || tolower(proc) == "y"){
    expr[,1] <- gsub("\\.\\d+$", "", expr[[1]])
    md = colnames(expr)[2:ncol(expr)]
    md %<>% as.data.frame()
    colnames(md) <- "id"
    
    n_groups <- as.numeric(readline(prompt = "Please enter the number of experimental conditions (excluding control): "))
    patterns <- vector(mode = "character", length = n_groups)
    groups <- vector(mode = "character", length = n_groups)
    ctrls <- vector(mode = "character", length = n_groups)
    
    for (i in 1:n_groups) {
      patterns[i] <- readline(prompt = paste("Please enter the unique identifier for experimental condition", i, ": "))
      groups[i] <- readline(prompt = paste("Please enter your desired name for experimental condition", i, ": "))
    }
    
    ctrls <- readline(prompt = paste("Please enter the name of your control: "))
    
    md %<>% mutate(condition = as.factor(
      case_when(
        grepl(patterns[1], id) ~ groups[1],
        TRUE ~ ctrls
      )
    ))
    
    for (i in 2:n_groups) {
      md %<>% mutate(condition = as.factor(
        case_when(
          grepl(patterns[i], id) ~ groups[i],
          TRUE ~ condition
        )
      ))
    }
    print("Done :>")
    return(md)
  }
  
  if(tolower(proc) == "no" || tolower(proc) == "n"){
    message("Please provide the dataframe with unique patterns for each experimental condition :)")
    cat("Alternatively, you can provide your own metadata file as the example above: \n")
  }
}

rundeseq <- function(expr, md){
  expr %<>% as.data.frame()
  suppressWarnings({
    gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    while(is.na(gene_col) == TRUE || gene_col > ncol(expr)){
      message("Error detected, please try again :)")
      gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    }
  })
  colnames(expr)[gene_col] <- "geneid"
  expr %<>% relocate(geneid)
  expr[, -1] <- lapply(expr[, -1], as.integer)
  suppressWarnings({
    sam_col <- as.numeric(readline(prompt = "Please enter the column number of your sample IDs in your metadata: "))
    while(is.na(sam_col) == TRUE || sam_col > ncol(md)){
      sam_col <- as.numeric(readline(prompt = "Please enter the column number of your sample IDs in your metadata: "))
    }
  })
  suppressWarnings({
    cdn_col <- as.numeric(readline(prompt = "Please enter the column number of your experimental condition in your metadata: "))
    while(is.na(cdn_col) == TRUE || cdn_col > ncol(md)){
      cdn_col <- as.numeric(readline(prompt = "Please enter the column number of your experimental condition in your metadata: "))
    }
  })
  colnames(md)[sam_col] <- "sampid"
  colnames(md)[cdn_col] <- "condition"
  expr %<>% dplyr::select(geneid, any_of(md$sampid))
  dds = DESeqDataSetFromMatrix(countData = expr,
                               colData = md,
                               design = ~condition, tidy = TRUE)
  ctrl <- readline(prompt = "Enter the control: ")
  dds$condition <- relevel(dds$condition, ref = ctrl)
  n_groups <- as.numeric(readline(prompt = "Enter the number of experimental conditions excluding control: "))
  
  prefil <- readline(prompt = "Do you want to pre-filter the count input? (yes/no): ")
  if (tolower(prefil) == "yes" || tolower(prefil) == "y"){
    count <- as.numeric(readline(prompt = "Enter your desired filter value. Default = 10: "))
    spl <- as.numeric(readline(prompt = "Please enter the number of the smallest biological replicates among your experimental conditions: "))
    keep <- rowSums(counts(dds) >= count) >= spl
    #keep <- rowSums(counts(dds)) >= count
    dds <- dds[keep,]
    
    rundds = DESeq(dds)
    
    p = as.numeric(readline(prompt = "Please enter your desired adjusted p-value between 0 and 1: "))
    
    if (n_groups == 1){
      print(head(results(rundds), tidy = TRUE))
      print(summary(results(rundds, alpha = p)))
      return(rundds) 
    }  
    else{
      nr <- as.numeric(resultsNames(rundds) %>% as.data.frame() %>% nrow())
      for (i in 2:nr){
        print(head(results(rundds, name = resultsNames(rundds)[i]), tidy = TRUE))
        print(summary(results(rundds, name = resultsNames(rundds)[i], alpha = p)))
      }
      return(rundds)
    }
  }
  else{
    rundds = DESeq(dds)
    
    p = as.numeric(readline(prompt = "Enter your desired adjusted p-value: "))
    
    if (n_groups == 1){
      print(head(results(rundds), tidy = TRUE))
      print(summary(results(rundds, alpha = p)))
      return(rundds) 
    }  
    else{
      nr <- as.numeric(resultsNames(rundds) %>% as.data.frame() %>% nrow())
      for (i in 2:nr){
        print(head(results(rundds, name = resultsNames(rundds)[i]), tidy = TRUE))
        print(summary(results(rundds, name = resultsNames(rundds)[i], alpha = p)))
      }
      return(rundds)
    }
  }
  cat("Analysis complete :>")
}

getresults <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }
  
  message("This function will produce the significant results in excel files based on your provided cutoffs")
  nr <- as.numeric(resultsNames(dds) %>% as.data.frame() %>% nrow())
  p = as.numeric(readline(prompt = "Please enter the desired adjusted p-value between 0 and 1: "))
  message("The log2 fold change(LFC) is usually between 0.6 to 1")
  lfc = as.numeric(readline(prompt = "Please enter the desired LFC cutoff: "))
  
  hugo <- readline(prompt = "Are the gene names provided in the expression data official HUGO gene symbols? (yes/no): ")
  if(hugo == "no" || hugo == "n"){
    src <- readline(prompt = "Are the gene names provided in the expression data Ensembl Gene IDs? (yes/no): ")
    if(src == "yes" || src == "y"){
      srcType = "ensembl_gene_id"
    }
    if(src == "no" || src == "n"){
      message("Please refer to the specific attributte name in attrlist")
      srcType <- readline(prompt = "Please enter the specific attribute name: ")
    }
  }
  
  if(hugo == "no"){message("Please note that genes where the official symbol are not found will be removed from the final result")
    res_hgnc <- as.data.frame(results(dds))
    res_hgnc$geneid <- sapply(strsplit(rownames(res_hgnc),split="\\+"), "[",1)
    if(srcType == "ensembl_gene_id"){res_hgnc$geneid <- gsub("\\.\\d+$", "", res_hgnc$geneid)}
    
    v <- res_hgnc$geneid
    cat("Please wait.....\n")
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
    
    ## Make sure there was at least one mapping
    if( nrow(ID) < 1 ) top_n( "No IDs mapped successfully" )
    
    ## Drop empty duds
    k <- which( ID[,2] == "" )
    if( length(k) > 0 ) ID <- ID[-k,]
    colnames(ID)[1] <- "geneid" }
  
  namelist <- function(x){
    names(x) <- c(resultsNames(dds)[-1])
    names(x) <- gsub("condition_", "", names(x))
    names(x) <- gsub("_", " ", names(x))
    return(x)
  }
  
  if (nr > 2){
    x <- as.numeric(nr-1)
    l1 <- vector("list", length = x)
    if(hugo == "no" || hugo == "n"){reml <- vector("list", length = x)}
    for(i in 1:x){
      ddsn <- i + 1
      cat("Now processing", resultsNames(dds)[ddsn], "\n")
      res = results(dds, name = resultsNames(dds)[ddsn], alpha = p)
      l1[[i]] <- res
      res %>% as.data.frame() -> resdf
      
      if(hugo == "yes"){resdf$hgnc_symbol <- rownames(resdf)}
      if(hugo == "no"){
        resdf$geneid <- sapply(strsplit(rownames(resdf),split="\\+"), "[",1)
        if(srcType == "ensembl_gene_id"){resdf$geneid <- gsub("\\.\\d+$", "", resdf$geneid)}
        rem_res <- resdf[!resdf$geneid %in% ID[,1], ]
        rem_res %<>% dplyr::select(-geneid)
        rem_res$geneid <- rownames(rem_res)
        resdf <- merge(resdf, ID, by = "geneid", all.x = TRUE)
        rownames(resdf) <- NULL
        reml[[i]] <- list(rem_res)
      }
      
      results = resdf[which(resdf$padj < p), ]
      resup = results[results$log2FoldChange > lfc, ]
      resdown = results[results$log2FoldChange < -abs(lfc), ]
      sigresults = rbind(resup, resdown)
      l1[[i]] <- list(res, resdf, sigresults)
      names(l1[[i]]) <- c("dds", "fullres", "sigres")
    }
    l1 <- namelist(l1)
    if(hugo == "no" || hugo == "n"){reml <- namelist(reml)}
    
    message("DESeq2 provides the function to compare the log fold change among the groups other than the control")
    com_res = readline(prompt = "Do you want the result comparisons among the query groups other than the control? (yes/no): ")
    if(tolower(com_res) == "yes" || tolower(com_res) == "y"){
     det_res = readline(prompt = "Among all your experimental control other than the control, do you want to compare all of them? (yes/no): ")
     
     if(det_res == "yes" || det_res == "y"){
       message("If so, please provide the metadata dataframe from the tidymd function ")
       md_name <- readline(prompt = "Please enter the name of the metadata dataframe from the tidymd function: ")
       tryCatch({
         md <- get(md_name)
       }, error = function(e){
         message("Error: ", e)
         message("Please provide a valid dataframe from the tidymd function: ")
         return(NULL)
       })
       condition_n <- as.numeric(readline(prompt = "Please enter the column number of the experimental condition: "))
       colnames(md)[condition_n] <- "condition"
       md$condition <- sapply(md$condition, as.factor)
       cdn <- md$condition
       query_ctrl <- readline(prompt = "Please enter the name of the control: ")
       cdn <- relevel(cdn, ref = query_ctrl)
       query <- levels(cdn)[-1]
       query_com <- combn(query, 2)
       qcn <- query_com %>% as.data.frame() %>% ncol()
       
       l2 <- vector("list", length = qcn)
       if(hugo == "no" || hugo == "n"){reml2 <- vector("list", length = qcn)}
       
       for (b in 1:qcn){
         query_group <- query_com[,b]
         cat("Now processing", query_group[1], "vs", query_group[2], "\n")
         res <- results(dds, contrast = c("condition", query_group), alpha = p)
         res %>% as.data.frame() -> resdf
         if(hugo == "yes"){resdf$hgnc_symbol <- rownames(resdf)}
         if(hugo == "no" || hugo =="n"){
           resdf$geneid <- sapply(strsplit(rownames(resdf),split="\\+"), "[",1)
           if(srcType == "ensembl_gene_id"){resdf$geneid <- gsub("\\.\\d+$", "", resdf$geneid)}
           rem_res <- resdf[!resdf$geneid %in% ID[,1], ]
           rem_res %<>% dplyr::select(-geneid)
           rem_res$geneid <- rownames(rem_res)
           resdf <- merge(resdf, ID, by = "geneid", all.x = TRUE)
           rownames(resdf) <- NULL
           reml2[[b]] <- list(rem_res)
         }
         results = resdf[which(resdf$padj < p), ]
         resup = results[results$log2FoldChange > lfc, ]
         resdown = results[results$log2FoldChange < -abs(lfc), ]
         sigresults = rbind(resup, resdown)
         l2[[b]] <- list(res, resdf, sigresults)
         names(l2[[b]]) <- c("dds", "fullres", "sigres")
        
         y <- query_group[1]
         z <- query_group[2]
         names(l2)[b] <- paste(y, "vs", z, sep = " ")
         if(hugo == "no" || hugo == "n"){names(reml2)[b] <- paste(y, "vs", z, sep = " ")}
       }
     }
     
     if(det_res == "no" || det_res == "n"){
       message("The following are your provided query groups: ")
       cdn <- as.character(unique(dds$condition))
       cdn_len <- as.numeric(length(cdn))
       
       for (n in 1:cdn_len){cat(paste(n, cdn[n], sep = "."), sep = "\n")}
       
       res_num <- as.numeric(readline(prompt = "How many results do you want to produce?: "))
       
       l2 <- vector("list", length = res_num)
       if(hugo == "no" || hugo == "n"){reml2 <- vector("list", length = res_num)}

       for (g in 1:res_num){
         cat("Now processing result", g, "\n")
         g1 <- readline(prompt = "Please enter the name of your first query group: ")
         g2 <- readline(prompt = "Please enter the name of your second query group: ")
         res = results(dds, contrast = c("condition", g1, g2), alpha = p)
         res %>% as.data.frame() -> resdf
         if(hugo == "yes"){resdf$hgnc_symbol <- rownames(resdf)}
         if(hugo == "no"){
           resdf$geneid <- sapply(strsplit(rownames(resdf),split="\\+"), "[",1)
           if(srcType == "ensembl_gene_id"){resdf$geneid <- gsub("\\.\\d+$", "", resdf$geneid)}
           rem_res <- resdf[!resdf$geneid %in% ID[,1], ]
           rem_res %<>% dplyr::select(-geneid)
           rem_res$geneid <- rownames(rem_res)
           resdf <- merge(resdf, ID, by = "geneid", all.x = TRUE)
           rownames(resdf) <- NULL
           reml2[[g]] <- list(rem_res)
         }
         results = resdf[which(resdf$padj < p), ]
         resup = results[results$log2FoldChange > lfc, ]
         resdown = results[results$log2FoldChange < -abs(lfc), ]
         sigresults = rbind(resup, resdown)
         l2[[g]] <- list(res, resdf, sigresults)
         names(l2[[g]]) <- c("dds", "fullres", "sigres")
         names(l2)[g] <- paste(g1, "vs", g2, sep = " ")
         if(hugo == "no" || hugo == "n"){names(reml2)[g] <- paste(g1, "vs", g2, sep = " ")}
       }
     }
     l1 <- c(l1, l2)
     if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- c(reml, reml2)}
    }
    save_df <- readline(prompt = "Do you want to save your results into individual excel files? (yes/no): ")
    if(tolower(save_df) == "yes" || tolower(save_df) == "y"){
      message("Which results do you want to save? All results(both), unfiltered results only(unfil), or significant results only(sig)?")
      save_res <- readline(prompt = "Please enter the results you want to save (both/unfil/sig): ")
      dir <- readline(prompt = "Do you want to create a new folder for the saved results? (yes/no): ")
      nlist <- as.numeric(length(l1))
  
      
      if(tolower(save_res) == "both"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            res2 <- as.data.frame(l1[[nl]][[3]])
            name1 <- paste("./", dir_n, "/", filen, "_fullres.xlsx", sep = "")
            name2 <- paste("./", dir_n, "/", filen, "_sigres.xlsx", sep = "")
            write_xlsx(res1, name1)
            write_xlsx(res2, name2)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            res2 <- as.data.frame(l1[[nl]][[3]])
            name1 <- paste("./", filen, "_fullres.xlsx", sep = "")
            name2 <- paste("./", filen, "_sigres.xlsx", sep = "")
            write_xlsx(res1, name1)
            write_xlsx(res2, name2)
          }
        }
      }
      if(tolower(save_res) == "unfil"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            name1 <- paste("./", dir_n, "/", filen, "_fullres.xlsx", sep = "")
            write_xlsx(res1, name1)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l2[[nl]][[2]])
            name1 <- paste("./", filen, "_fullres.xlsx", sep = "")
            write_xlsx(res1, name1)
          }
        }
      }
      if(tolower(save_res) == "sig"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res2 <- as.data.frame(l1[[nl]][[3]])
            name2 <- paste("./", dir_n, "/", filen, "_sigres.xlsx", sep = "")
            write_xlsx(res2, name2)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res2 <- as.data.frame(l1[[nl]][[3]])
            name2 <- paste("./", filen, "_sigres.xlsx", sep = "")
            write_xlsx(res2, name2)
          }
        }
      }
      if(tolower(hugo) == "no" || tolower(hugo) == "n"){
        del_genes <- readline(prompt = "Do you want the dataframes of the deleted genes where the official gene symbol was not found? (yes/no): ")
        if(tolower(del_genes) == "yes" || tolower(del_genes) == "y"){
          if(tolower(dir) == "yes" || tolower(dir) == "y"){
            for(dl in seq_along(reml)){
              dfilen <- names(reml)[dl]
              cat("Now saving", dfilen, "\n")
              dres <- as.data.frame(reml[[dl]])
              name1 <- paste("./", dir_n, "/", dfilen, "_removedgenes.xlsx", sep = "")
              write_xlsx(dres, name1)
            }
          }
          if(tolower(dir) == "no" || tolower(dir) == "n"){
            for(dl in seq_along(reml)){
              dfilen <- names(reml)[dl]
              cat("Now saving", dfilen, "\n")
              dres <- as.data.frame(reml[[dl]])
              name1 <- paste("./", dfilen, "_removedgenes.xlsx", sep = "")
              write_xlsx(dres, name1)
            }
          }
        }
      }
    }
    return(l1)
  }
  else{
    res = results(dds, alpha = p)
    res %>% as.data.frame() -> resdf
    if(hugo == "no"){
      resdf$geneid <- sapply(strsplit(rownames(resdf),split="\\+"), "[",1)
      if(srcType == "ensembl_gene_id"){resdf$geneid <- gsub("\\.\\d+$", "", resdf$geneid)}
      rem_res <- resdf[!resdf$geneid %in% ID[,1], ]
      rem_res %<>% dplyr::select(-geneid)
      rem_res$geneid <- rownames(rem_res)
      resdf <- merge(resdf, ID, by = "geneid", all.x = TRUE)
      rownames(resdf) <- NULL
    }
    if(hugo == "yes"){resdf$hgnc_symbol <- rownames(resdf)}
    results = resdf[which(resdf$padj < p), ]
    resup = results[results$log2FoldChange > lfc, ]
    resdown = results[results$log2FoldChange < -abs(lfc), ]
    sigresults = rbind(resup, resdown)
    
    filen <- readline(prompt = "Please provide the name of your result: ")
    name1 <- paste("./", filen, "_fullres.xlsx", sep = "")
    name2 <- paste("./", filen, "_sigres.xlsx", sep = "")
    name3 <- paste("./", filen, "_deletedgenes.xlsx", sep = "")
    write_xlsx(resdf, name1)
    write_xlsx(sigresults, name2)
    if(hugo == "no"){write_xlsx(rem_res, name3)}
    
    return(res)
  }
  cat("Done :>")
}

shrinklfc <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }
  
  readmd <- function(){
    message("If so, please provide the metadata dataframe from the tidymd function ")
    md_name <- readline(prompt = "Please enter the name of the metadata dataframe from the tidymd function: ")
    tryCatch({
      md <- get(md_name)
      
    }, error = function(e){
      message("Error: ", e)
      message("Please provide a valid dataframe from the tidymd function: ")
      return(NULL)
    })
    return(md)
  }
  namelist <- function(x){
    names(x) <- c(resultsNames(dds)[-1])
    names(x) <- gsub("condition_", "", names(x))
    names(x) <- gsub("_", " ", names(x))
    return(x)
  }
  
  
  p <- as.numeric(readline(prompt = "Please enter your desired adjusted p-value: "))
  lfc <- as.numeric(readline(prompt = "Please enter your desired LFC cutoff: "))
  x = as.numeric(length(resultsNames(dds)))
  
  hugo <- readline(prompt = "Are the gene names provided in the expression data official HUGO gene symbols? (yes/no): ")
  if(hugo == "no" || hugo == "n"){
    src <- readline(prompt = "Are the gene names provided in the expression data Ensembl IDs? (yes/no): ")
    if(src == "yes" || src == "y"){
      srcType = "ensembl_gene_id"
    }
    if(src == "no" || src == "n"){
      message("Please refer to the specific attributte name in attrlist")
      srcType <- readline(prompt = "Please enter the specific attribute name: ")
    }
  }
  
  if(hugo == "no"){message("Please note that genes where the official symbol are not found will be removed from the final result")
    res_hgnc <- as.data.frame(results(dds))
    res_hgnc$geneid <- sapply(strsplit(rownames(res_hgnc),split="\\+"), "[",1)
    if(srcType == "ensembl_gene_id"){res_hgnc$geneid <- gsub("\\.\\d+$", "", res_hgnc$geneid)}
    
    v <- res_hgnc$geneid
    cat("Please wait.....\n")
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
    
    ## Make sure there was at least one mapping
    if( nrow(ID) < 1 ) top_n( "No IDs mapped successfully" )
    
    ## Drop empty duds
    k <- which( ID[,2] == "" )
    if( length(k) > 0 ) ID <- ID[-k,]
    colnames(ID)[1] <- "geneid" }
  
  if (x > 2){
    ans1 <- readline(prompt = "Do you want to shrink all the results? (yes/no): ")
    message("Please note that the default estimator (apeglm) will only shrink results of your experimental conditions against the control")
    Sys.sleep(0.5)
    message("If you want to shrink the results comparing the exprimental conditions other than the control, \nplease change the shrinkage esimator to ashr/normal")
    message("However, do note that according to the author, apeglm and ashr tend to show less bias than normal")
    Sys.sleep(0.5)
    message("For more information about LFC shrinkage, please refer to the DESeq2 vignette")
    message("http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html")
    Sys.sleep(0.5)
    ans2 <- readline(prompt = "Do you want to change the shrinkage estimator? Default = apeglm (yes/no): ")
    
    if (tolower(ans1) == "yes" || tolower(ans1) == "y"){
      y <- as.numeric(x - 1)
      
      if(tolower(ans2) == "yes" || tolower(ans2) == "yes"){
        ans3 <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
        com_res <- readline(prompt = "Do you want to shrink your results among your experimental conditions, excluding the control? (yes/no): ")
        
        if(tolower(com_res) == "yes" || tolower(com_res) == "y"){
          md <- readmd()
          condition_n <- as.numeric(readline(prompt = "Please enter the column number of the experimental condition: "))
          colnames(md)[condition_n] <- "condition"
          md$condition <- sapply(md$condition, as.factor)

          cdn <- md$condition
          query_ctrl <- readline(prompt = "Please enter the name of the control: ")
          cdn <- relevel(cdn, ref = query_ctrl)
          query <- levels(cdn)[-1]
          query_com <- combn(query, 2)
          qcn <- as.numeric(query_com %>% as.data.frame() %>% ncol())
          
          flength <- y + qcn
          l1 <- vector("list", length = y)
          if (hugo == "no" || hugo == "n"){reml <- vector("list", length = y)}
          for(i in 1:y){
            resn <- as.numeric(i + 1)
            resname <- resultsNames(dds)
            resname <- gsub("condition_", "", resname)
            resname <- gsub("_", " ", resname)
            cat("Now shrinking result ", i, ": ", resname[resn], "\n", sep = "")
            slr <- lfcShrink(dds, coef = resn, type = ans3, res = results(dds, name = resultsNames(dds)[resn], alpha = p))
            slr %>% as.data.frame() -> slrdf
            if(hugo == "yes"){
              slrdf$hgnc_symbol <- rownames(slrdf)
              rownames(slrdf) <- NULL
              }
            if(hugo == "no"){
              slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
              if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
              rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
              rem_res %<>% dplyr::select(-geneid)
              rem_res$geneid <- rownames(rem_res)
              slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
              rownames(slrdf) <- NULL
              reml[[i]] <- list(rem_res)
            }
            
            sig_slr = slrdf[which(slrdf$padj < p), ]
            slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
            slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
            sigresults = rbind(slr_up, slr_down)
            l1[[i]] <- list(slr, slrdf, sigresults)
            names(l1[[i]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
          }
          l1 <- namelist(l1)
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- namelist(reml)}
          
          l2 <- vector("list", length = qcn)
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml2 <- vector("list", length = qcn)}
          for(g in 1:qcn){
            resn2 <- as.numeric(y + g)
            query_group <- query_com[,g]
            n1 <- query_group[1]
            n2 <- query_group[2]
            cat("Now shrinking result ", resn2, ": ", n1, " vs ", n2, "\n", sep = "")
            slr2 <- lfcShrink(dds, contrast = c("condition", query_group), type = ans3, alpha = p)
            slr2 %>% as.data.frame() -> slrdf
            if(hugo == "yes"){
              slrdf$hgnc_symbol <- rownames(slrdf)
              rownames(slrdf) <- NULL
            }
            if(hugo == "no"){
              slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
              if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
              rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
              rem_res %<>% dplyr::select(-geneid)
              rem_res$geneid <- rownames(rem_res)
              slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
              rownames(slrdf) <- NULL
              reml2[[g]] <- list(rem_res)
            }
            sig_slr = slrdf[which(slrdf$padj < p), ]
            slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
            slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
            sigresults = rbind(slr_up, slr_down)
            l2[[g]] <- list(slr2, slrdf, sigresults)
            names(l2[[g]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
            names(l2)[[g]] <- paste(n1, "vs", n2, sep = " ")
            if(tolower(hugo) == "no" || tolower(hugo) == "n"){names(reml2)[[g]] <- paste(n1, "vs", n2, sep = " ")}
          }
          l1 <- c(l1,l2)
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- c(reml, reml2)}
        }
        
        if(tolower(com_res) == "no" || tolower(com_res) == "n"){
          l1 <- vector("list", length = y)
          for(i in 1:y){
            resn <- as.numeric(i + 1)
            resname <- resultsNames(dds)
            resname <- gsub("condition_", "", resname)
            resname <- gsub("_", " ", resname)
            cat("Now shrinking result ", i, ": ", resname[resn], "\n", sep = "")
            slr <- lfcShrink(dds, coef = resn, type = ans3, res = results(dds, name = resultsNames(dds)[resn], alpha = p))
            slr %>% as.data.frame() -> slrdf
            if(hugo == "yes"){
              slrdf$hgnc_symbol <- rownames(slrdf)
              rownames(slrdf) <- NULL
            }
            if(hugo == "no"){
              slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
              if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
              rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
              rem_res %<>% dplyr::select(-geneid)
              rem_res$geneid <- rownames(rem_res)
              slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
              rownames(slrdf) <- NULL
              reml[[i]] <- list(rem_res)
            }
            
            sig_slr = slrdf[which(slrdf$padj < p), ]
            slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
            slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
            sigresults = rbind(slr_up, slr_down)
            l1[[i]] <- list(slr, slrdf, sigresults)
            names(l1[[i]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
          }
          l1 <- namelist(l1)
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- namelist(reml)}
        }
      }
      
      if (tolower(ans2) == "no" || tolower(ans2) == "n"){
        l1 <- vector("list", length = y)
        if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- vector("list", length = y)}
        for(i in 1:y){
          resn <- as.numeric(i + 1)
          resname <- resultsNames(dds)
          resname <- gsub("condition_", "", resname)
          resname <- gsub("_", " ", resname)
          cat("Now shrinking result ", i, ": ", resname[resn], "\n", sep = "")
          slr <- lfcShrink(dds, coef = resn, type = "apeglm", res = results(dds, name = resultsNames(dds)[resn], alpha = p))
          slr %>% as.data.frame() -> slrdf
          if(hugo == "yes"){
            slrdf$hgnc_symbol <- rownames(slrdf)
            rownames(slrdf) <- NULL
          }
          if(hugo == "no"){
            slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
            if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
            rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
            rem_res %<>% dplyr::select(-geneid)
            rem_res$geneid <- rownames(rem_res)
            slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
            rownames(slrdf) <- NULL
            reml[[i]] <- list(rem_res)
          }
          
          sig_slr = slrdf[which(slrdf$padj < p), ]
          slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
          slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
          sigresults = rbind(slr_up, slr_down)
          l1[[i]] <- list(slr, slrdf, sigresults)
          names(l1[[i]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
        }
        l1 <- namelist(l1)
        if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- namelist(reml)}
      }
    }
    if(tolower(ans1) == "no" || tolower(ans1) == "n"){
      res_num <- as.numeric(readline(prompt = "How many results do you want to shrink?: "))
      l1 <- vector("list", length = res_num)
      if(tolower(hugo) == "no" || tolower(hugo) == "n"){reml <- vector("list", length = res_num)}
      if(res_num == 1){
        if (tolower(ans2) == "no" || tolower(ans2) == "n"){
          cat("The following are your listed conditions (starting from number 2) :) \n")
          dds_con <- resultsNames(dds)
          dds_con <- gsub("condition_", "", dds_con)
          dds_con <- gsub("_", " ", dds_con)
          for (dc in 1:length(resultsNames(dds))) {cat(dc, ". ", dds_con[dc], "\n", sep = "")}
          c <- as.numeric(readline(prompt = "Enter the number of your desired query condition: "))
          slr <- lfcShrink(dds, coef = c, type = "apeglm", res = results(dds, name = resultsNames(dds)[c], alpha = p))
          slr %>% as.data.frame() -> slrdf
          if(hugo == "yes"){
            slrdf$hgnc_symbol <- rownames(slrdf)
            rownames(slrdf) <- NULL
          }
          if(hugo == "no"){
            slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
            if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
            rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
            rem_res %<>% dplyr::select(-geneid)
            rem_res$geneid <- rownames(rem_res)
            slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
            rownames(slrdf) <- NULL
            reml[[1]] <- list(rem_res)
          }
          
          sig_slr = slrdf[which(slrdf$padj < p), ]
          slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
          slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
          sigresults = rbind(slr_up, slr_down)
          l1[[1]] <- list(slr, slrdf, sigresults)
          names(l1[[1]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
          lname <- resultsNames(dds)[c]
          lname <- gsub("condition_", "", lname)
          lname <- gsub("_", " ", lname)
          names(l1) <- lname
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){names(reml) <- lname}
        }
        if (tolower(ans2) == 'yes' || tolower(ans2) == "y") {
          ans4 <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
          cond <- as.character(unique(dds$condition))
          message("The following are your listed experimental conditions: ")
          for (v in 1:length(cond)){cat(v, ".", cond[v], "\n", sep = "")}
          c1 <- readline(prompt = "Please enter your first experimental condition: ")
          c2 <- readline(prompt = "Please enter your second experimental condtion: ")
          slr <- lfcShrink(dds, type = ans4, contrast = c("condition", c1, c2), alpha = p)
          slr %>% as.data.frame() -> slrdf
          if(hugo == "yes"){
            slrdf$hgnc_symbol <- rownames(slrdf)
            rownames(slrdf) <- NULL
          }
          if(hugo == "no"){
            slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
            if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
            rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
            rem_res %<>% dplyr::select(-geneid)
            rem_res$geneid <- rownames(rem_res)
            slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
            rownames(slrdf) <- NULL
            reml[[1]] <- list(rem_res)
          }
          
          sig_slr = slrdf[which(slrdf$padj < p), ]
          slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
          slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
          sigresults = rbind(slr_up, slr_down)
          l1[[1]] <- list(slr, slrdf, sigresults)
          names(l1[[1]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
          lname <- resultsNames(dds)[c]
          lname <- gsub("condition_", "", lname)
          lname <- gsub("_", " ", lname)
          names(l1) <- lname
          if(tolower(hugo) == "no" || tolower(hugo) == "n"){names(reml) <- lname}
        }
      }
      if(res_num >= 2){
        if(tolower(ans2) == "no" || tolower(ans2) == "n"){
          cat("The following are your listed conditions (starting from number 2) :) \n")
          dds_con <- resultsNames(dds)
          dds_con <- gsub("condition_", "", dds_con)
          dds_con <- gsub("_", " ", dds_con)
          for (dc in 1:length(resultsNames(dds))) {cat(dc, ". ", dds_con[dc], "\n", sep = "")}
          for(k in 1:res_num){
            cat("Now processing result", k, "\n")
            c <- as.numeric(readline(prompt = "Please enter the number of your desired query condition: "))
            plotname <- resultsNames(dds)[c]
            slr <- lfcShrink(dds, coef = c, type = "apeglm", res = results(dds, name = resultsNames(dds)[c], alpha = p))
            slr %>% as.data.frame() -> slrdf
            if(hugo == "yes"){
              slrdf$hgnc_symbol <- rownames(slrdf)
              rownames(slrdf) <- NULL
            }
            if(hugo == "no"){
              slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
              if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
              rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
              rem_res %<>% dplyr::select(-geneid)
              rem_res$geneid <- rownames(rem_res)
              slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
              rownames(slrdf) <- NULL
              reml[[k]] <- list(rem_res)
            }
            
            sig_slr = slrdf[which(slrdf$padj < p), ]
            slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
            slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
            sigresults = rbind(slr_up, slr_down)
            l1[[k]] <- list(slr, slrdf, sigresults)
            names(l1[[k]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
            plotname <- gsub("condition_", "", plotname)
            plotname <- gsub("_", " ", plotname)
            names(l1)[k] <- plotname
            if(tolower(hugo) == "no" || tolower(hugo) == "n"){names(reml)[k] <- plotname}
          }
        }
        if(tolower(ans2) == "yes" || tolower(ans2) == "y"){
          ans4 <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
          cond <- as.character(unique(dds$condition))
          l1 <- vector("list", length = res_num)
          message("The following are your listed experimental conditions: ")
          for (v in 1:length(cond)){cat(v, ".", cond[v], "\n", sep = "")}
          for(j in 1:res_num){
            cat("Now processing result", j, "\n")
            c1 <- readline(prompt = "Please enter your first experimental condition: ")
            c2 <- readline(prompt = "Please enter your second experimental condtion: ")
            slr <- lfcShrink(dds, type = ans4, contrast = c("condition", c1, c2), alpha = p)
            slr %>% as.data.frame() -> slrdf
            if(hugo == "yes"){
              slrdf$hgnc_symbol <- rownames(slrdf)
              rownames(slrdf) <- NULL
            }
            if(hugo == "no"){
              slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
              if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
              rem_res <- slrdf[!slrdf$geneid %in% ID[,1], ]
              rem_res %<>% dplyr::select(-geneid)
              rem_res$geneid <- rownames(rem_res)
              slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
              rownames(slrdf) <- NULL
              reml[[j]] <- list(rem_res)
            }
            
            sig_slr = slrdf[which(slrdf$padj < p), ]
            slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
            slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
            sigresults = rbind(slr_up, slr_down)
            l1[[j]] <- list(slr, slrdf, sigresults)
            names(l1[[j]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
            names(l1)[j] <- paste(c1, "vs", c2, sep = " ")
            if(tolower(hugo) == "no" || tolower(hugo) == "n"){names(reml)[j] <- paste(c1, "vs", c2, sep = " ")}
          }
        }
      }
    }
    save_df <- readline(prompt = "Do you want to save your results into individual excel files? (yes/no): ")
    if(tolower(save_df) == "yes" || tolower(save_df) == "y"){
      message("Which results do you want to save? All results(both), unfiltered results only(unfil), or significant results only(sig)?")
      save_res <- readline(prompt = "Please enter the results you want to save (both/unfil/sig): ")
      dir <- readline(prompt = "Do you want to create a new folder for the saved results? (yes/no): ")
      nlist <- as.numeric(length(l1))
      
      
      if(tolower(save_res) == "both"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            res2 <- as.data.frame(l1[[nl]][[3]])
            name1 <- paste("./", dir_n, "/", filen, "_shrinked_fullres.xlsx", sep = "")
            name2 <- paste("./", dir_n, "/", filen, "_shrinked_sigres.xlsx", sep = "")
            write_xlsx(res1, name1)
            write_xlsx(res2, name2)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            res2 <- as.data.frame(l1[[nl]][[3]])
            name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
            name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
            write_xlsx(res1, name1)
            write_xlsx(res2, name2)
          }
        }
      }
      if(tolower(save_res) == "unfil"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l1[[nl]][[2]])
            name1 <- paste("./", dir_n, "/", filen, "_shrinked_fullres.xlsx", sep = "")
            write_xlsx(res1, name1)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res1 <- as.data.frame(l2[[nl]][[2]])
            name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
            write_xlsx(res1, name1)
          }
        }
      }
      if(tolower(save_res) == "sig"){
        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          message("Please make sure the name of your created folder does not exist in your current directory")
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          dir.create(dir_n)
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res2 <- as.data.frame(l1[[nl]][[3]])
            name2 <- paste("./", dir_n, "/", filen, "_shrinked_sigres.xlsx", sep = "")
            write_xlsx(res2, name2)
          }
        }
        if(tolower(dir) == "no" || tolower(dir) == "n"){
          for (nl in seq_along(l1)){
            filen <- names(l1)[nl]
            cat("Now saving", filen, "\n")
            res2 <- as.data.frame(l1[[nl]][[3]])
            name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
            write_xlsx(res2, name2)
          }
        }
      }
      if(tolower(hugo) == "no" || tolower(hugo) == "n"){
        del_genes <- readline(prompt = "Do you want the dataframes of the deleted genes where the official gene symbol was not found? (yes/no): ")
        if(tolower(del_genes) == "yes" || tolower(del_genes) == "y"){
          if(tolower(dir) == "yes" || tolower(dir) == "y"){
            for(dl in seq_along(reml)){
              dfilen <- names(reml)[dl]
              cat("Now saving", dfilen, "\n")
              dres <- as.data.frame(reml[[dl]])
              name1 <- paste("./", dir_n, "/", dfilen, "_shrinked_removedgenes.xlsx", sep = "")
              write_xlsx(dres, name1)
            }
          }
          if(tolower(dir) == "no" || tolower(dir) == "n"){
            for(dl in seq_along(reml)){
              dfilen <- names(reml)[dl]
              cat("Now saving", dfilen, "\n")
              dres <- as.data.frame(reml[[dl]])
              name1 <- paste("./", dfilen, "_shrinked_removedgenes.xlsx", sep = "")
              write_xlsx(dres, name1)
            }
          }
        }
      }
    }
  }
  if (x == 2){
    ans5 <- readline(prompt = "Do you want to change the shrinkage estimator? Default = apeglm (yes/no): ")
    l1 <- vector("list", length = 1)
    if (tolower(ans5) == "no" || tolower(ans5) == "n"){
      slr <- lfcShrink(dds, coef = 2, type = "apeglm", res = results(dds, name = resultsNames(dds)[2], alpha = p))
      slr %>% as.data.frame() -> slrdf
      if(hugo == "yes"){
        slrdf$hgnc_symbol <- rownames(slrdf)
        rownames(slrdf) <- NULL
      }
      if(hugo == "no"){
        slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
        if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
        slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
        rownames(slrdf) <- NULL
      }
      
      sig_slr = slrdf[which(slrdf$padj < p), ]
      slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
      slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
      sigresults = rbind(slr_up, slr_down)
      l1[[1]] <- list(slr, slrdf, sigresults)
      names(l1[[1]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
      lname <- resultsNames(dds)[2]
      lname <- gsub("condition_", "", lname)
      lname <- gsub("_", " ", lname)
      names(l1) <- lname
      res1 <- as.data.frame(l1[[1]][[2]])
      res2 <- as.data.frame(l1[[1]][[3]])
      filen <- readline(prompt = "Please provide a name for your files: ")
      name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
      name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
      write_xlsx(res1, name1)
      write_xlsx(res2, name2)
    }
    if (tolower(ans5) == "yes" || tolower(ans5) == "y") {
      ans6 <- readline(prompt = "Please enter your desired shrinkage estimator(apeglm/normal/ashr): ")
      slr <- lfcShrink(dds, coef = 2, type = ans6, res = results(dds, name = resultsNames(dds)[2], alpha = p))
      slr %>% as.data.frame() -> slrdf
      if(hugo == "yes"){
        slrdf$hgnc_symbol <- rownames(slrdf)
        rownames(slrdf) <- NULL
      }
      if(hugo == "no"){
        slrdf$geneid <- sapply(strsplit(rownames(slrdf),split="\\+"), "[",1)
        if(srcType == "ensembl_gene_id"){slrdf$geneid <- gsub("\\.\\d+$", "", slrdf$geneid)}
        slrdf <- merge(slrdf, ID, by = "geneid", all.x = TRUE)
        rownames(slrdf) <- NULL
      }
      
      sig_slr = slrdf[which(slrdf$padj < p), ]
      slr_up = sig_slr[sig_slr$log2FoldChange > lfc,]
      slr_down = sig_slr[sig_slr$log2FoldChange < -abs(lfc), ]
      sigresults = rbind(slr_up, slr_down)
      l1[[1]] <- list(slr, slrdf, sigresults)
      names(l1[[1]]) <- c("lfcdds", "_shrinked_fullres", "shrinked_sigres")
      lname <- resultsNames(dds)[2]
      lname <- gsub("condition_", "", lname)
      lname <- gsub("_", " ", lname)
      names(l1) <- lname
      res1 <- as.data.frame(l1[[1]][[2]])
      res2 <- as.data.frame(l1[[1]][[3]])
      filen <- readline(prompt = "Please provide a name for your files: ")
      name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
      name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
      write_xlsx(res1, name1)
      write_xlsx(res2, name2)
    }
  }
  cat("Done :>")
  return(l1)
}

maplot <- function(result){
  if(!inherits(result, c("DESeqResults", "list"))){
    message("Please provide the DESeqResults from the getresults/shrinklfc function")
    return(NULL)
  }
  if(class(result) == "list"){
    for (m in seq_along(result)){
      cs <- class(result[[m]][[1]]) == "DESeqResults"
      if (cs != "TRUE"){
        message("Please provide the DESeqResults from the getresults/shrinklfc function")
        return(NULL)
      }
    }
  }
  
  def_map <- function(x){
    for (l in 1:length(x)){
      i = x[l]
      cat("Now processing", names(result)[i], "\n")
      plot <- plotMA(result[[i]][[1]], ylim = c(-10, 10), 
                     xlab = "Mean of Normalised Counts", 
                     ylab = "Log2 Fold Change", 
                     main = names(result)[i])
    }
  }
  def_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    for (l in 1:length(x)){
      i = x[l]
      cat("Now processing", names(result)[i], "\n")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        name <- paste0(dir_n, "/", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      if(tolower(dir) == "no" || tolower(dir) == "n"){
        name <- paste0("./", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      plot <- plotMA(result[[i]][[1]], ylim = c(-10, 10), 
                     xlab = "Mean of Normalised Counts", 
                     ylab = "Log2 Fold Change", 
                     main = names(result)[i])
      dev.off()
    }
  } 
  
  def_title_map <- function(x){
    for (l in 1:length(x)){
      i = x[l]
      cat("Now processing plot", i, "\n")
      cat(names(result)[i], "\n")
      title <- readline(prompt = "Please provide the title for the current plot: ")
      plot <- plotMA(result[[i]][[1]], ylim = c(-10, 10), 
                     xlab = "Mean of Normalised Counts", 
                     ylab = "Log2 Fold Change", 
                     main = title)
    }
  }
  def_title_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    for (l in 1:length(x)){
      i = x[l]
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        name <- paste0(dir_n, "/", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      if(tolower(dir) == "no" || tolower(dir) == "n"){
        name <- paste0("./", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      cat("Now processing plot", i, "\n")
      cat(names(result)[i], "\n")
      title <- readline(prompt = "Please provide the title for the current plot: ")
      plot <- plotMA(result[[i]][[1]], ylim = c(-10, 10), 
                     xlab = "Mean of Normalised Counts", 
                     ylab = "Log2 Fold Change", 
                     main = title)
      dev.off()
    }
  }
  
  ndef_map <- function(x){
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10: "))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    for (l in 1:length(x)){
      i = x[l]
      cat("Now processing", names(result)[i], "\n")
      plot <- plotMA(result[[i]][[1]], ylim = c(-abs(lowerlimit), upperlimit), 
                     xlab = labx, 
                     ylab = laby, 
                     main = names(result)[i])
    }
  }
  ndef_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10:"))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    for (l in 1:length(x)){
      i = x[l]
      cat("Now processing", names(result)[i], "\n")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        name <- paste0(dir_n, "/", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      if(tolower(dir) == "no" || tolower(dir) == "n"){
        name <- paste0("./", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      plot <- plotMA(result[[i]][[1]], ylim = c(-abs(lowerlimit), upperlimit), 
                     xlab = labx, 
                     ylab = laby, 
                     main = names(result)[i])
      dev.off()
    }
  }
  
  ndef_title_map <- function(x){
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10: "))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    for (i in 1:length(x)){
      i = x[l]
      cat("Now processing plot", i, "\n")
      cat(names(result)[i], "\n")
      title <- readline(prompt = "Please provide the title for the current plot: ")
      plot <- plotMA(result[[i]][[1]], ylim = c(-abs(lowerlimit), upperlimit), 
                     xlab = labx, 
                     ylab = laby, 
                     main = title)
    }
  }
  ndef_title_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10: "))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    for (i in 1:length(x)){
      i = x[l]
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        name <- paste0(dir_n, "/", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      if(tolower(dir) == "no" || tolower(dir) == "n"){
        name <- paste0("./", names(result)[i], ".png", sep = "")
        png(name, res = 300, units = "in", width = 8, height = 6)
      }
      cat("Now processing plot", i, "\n")
      cat(names(result)[i], "\n")
      title <- readline(prompt = "Please provide the title for the current plot: ")
      plot <- plotMA(result[[i]][[1]], ylim = c(-abs(lowerlimit), upperlimit), 
                     xlab = labx, 
                     ylab = laby, 
                     main = title)
      dev.off()
    }
  }
  
  dsr_def_title_map <- function(x){
    title <- readline(prompt = "Please provide the title for the current plot: ")
    plot <- plotMA(result, ylim = c(-10, 10), 
                   xlab = "Mean of Normalised Counts", 
                   ylab = "Log2 Fold Change", 
                   main = title)
  }
  dsr_def_title_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    title <- readline(prompt = "Please provide the title for the current plot: ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      name <- paste0(dir_n, "/", title, ".png", sep = "")
      png(name, res = 300, units = "in", width = 8, height = 6)
    }
    if(tolower(dir) == "no" || tolower(dir) == "n"){
      name <- paste0("./", title, ".png", sep = "")
      png(name, res = 300, units = "in", width = 8, height = 6)
    }
    plot <- plotMA(result, ylim = c(-10, 10), 
                   xlab = "Mean of Normalised Counts", 
                   ylab = "Log2 Fold Change", 
                   main = title)
    dev.off()
  }
  
  dsr_ndef_title_map <- function(x){
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10: "))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    title <- readline(prompt = "Please provide the title for the current plot: ")
    plot <- plotMA(result, ylim = c(-abs(lowerlimit), upperlimit), 
                   xlab = labx, 
                   ylab = laby, 
                   main = title)
  }
  dsr_ndef_title_save_map <- function(x){
    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      message("Please make sure the name of your created folder does not exist in your current directory")
      dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
      dir.create(dir_n)
    }
    title <- readline(prompt = "Please provide the title for the current plot: ")
    lowerlimit = as.numeric(readline(prompt = "Please enter the value of the lower limit. Default = 10: "))
    upperlimit = as.numeric(readline(prompt = "Please enter the value of the upper limit. Default = 10: "))
    labx = readline(prompt = "Please enter your desired x-axis label: ")
    laby = readline(prompt = "Please enter your desired y-axis label: ")
    
    if(tolower(dir) == "yes" || tolower(dir) == "y"){
      name <- paste0(dir_n, "/", title, ".png", sep = "")
      png(name, res = 300, units = "in", width = 8, height = 6)
    }
    if(tolower(dir) == "no" || tolower(dir) == "n"){
      name <- paste0("./", title, ".png", sep = "")
      png(name, res = 300, units = "in", width = 8, height = 6)
    }
    plot <- plotMA(result, ylim = c(-abs(lowerlimit), upperlimit), 
                   xlab = labx, 
                   ylab = laby, 
                   main = title)
    dev.off()
  }
  
  if(class(result) == "list"){
    all_res <- readline(prompt = "Do you want to plot all your results? (yes/no): ")
    
    if(tolower(all_res) == "yes" || tolower(all_res) == "y"){
      nplot <- as.numeric(seq(1, length(result)))
      def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
      ans_title <- readline(prompt = "Do you want to supply your own plot titles? (yes/no): ")
      save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
      if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){def_title_save_map(nplot)}
      if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){def_title_map(nplot)}
      if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){def_save_map(nplot)}
      if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "no"){def_map(nplot)}
      if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_map(nplot)}
      if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){ndef_title_map(nplot)}
      if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){ndef_save_map(nplot)}
      if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "no"){ndef_map(nplot)}
    }
    
    if(tolower(all_res) == "no" || tolower(all_res) == "n"){
      message("The following are your provided results: ")
      for (np in 1:length(result)){cat(np, ". ", names(result)[np], "\n", sep = "")}
      num_plot <- as.numeric(readline(prompt = "Please state the number of plots required: "))
      if(num_plot >= 2){
        message("Please specify the result number for the plots in a comma separated format WITH ONE SPACES: e.g. 1, 2, 4")
        splot <- readline(prompt = "Please enter the result number of the plots required: ")
        splotn <- as.numeric(strsplit(splot, ", ")[[1]])
        def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
        ans_title <- readline(prompt = "Do you want to supply your own plot titles? (yes/no): ")
        save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
        if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){def_title_save_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){def_title_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){def_save_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "no"){def_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){ndef_title_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){ndef_save_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "no"){ndef_map(splotn)}
      }
      if(num_plot == 1){
        splotn <- as.numeric(readline(prompt = "Please specify the result number for the plot required: "))
        def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
        ans_title <- readline(prompt = "Do you want to supply your own plot titles? (yes/no): ")
        save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
        if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){def_title_save_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){def_title_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){def_save_map(splotn)}
        if(tolower(def_par) == "yes" && tolower(ans_title) == "no" && tolower(save_res) == "no"){def_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "yes" && tolower(save_res) == "no"){ndef_title_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "yes"){ndef_save_map(splotn)}
        if(tolower(def_par) == "no" && tolower(ans_title) == "no" && tolower(save_res) == "no"){ndef_map(splotn)}
      }
    }
  }
  if(class(result) == "DESeqResults"){
    def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
    save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
    if(tolower(def_par) == "yes" && tolower(save_res) == "yes"){dsr_def_title_save_map(result)}
    if(tolower(def_par) == "yes" && tolower(save_res) == "no"){dsr_def_title_map(result)}
    if(tolower(def_par) == "no" && tolower(save_res) == "yes"){dsr_ndef_title_save_map(result)}
    if(tolower(def_par) == "no" && tolower(save_res) == "no"){dsr_ndef_title_map(result)}
  }
}

volcanoplot <- function(result){
  if(!inherits(result, c("DESeqResults", "list"))){
    message("Please provide the DESeqResults from the getresults/shrinklfc function")
    return(NULL)
  }
  if(class(result) == "list"){
    for (m in seq_along(result)){
      cs <- class(result[[m]][[1]]) == "DESeqResults"
      if (cs != "TRUE"){
        message("Please provide the DESeqResults from either the getresults or shrinklfc function")
        return(NULL)
      }
    }
  }
  slab_func <- function(){
    message("Do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
    input = readline(prompt = "Enter D or V: ")
    if(tolower(input) == "d"){
      message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
      genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
      slabel <- strsplit(genes, ", ")[[1]]
    }
    else if(tolower(input) == 'v'){
      message("Please provide a dataframe or a vector containing the gene names :)")
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        slabel <- get(genes)
        
      }, error = function(e){
        message("Error: ", e)
        message("Please provide a valid dataframe or a vector of gene names: ")
        return(NULL)
      })
      
      cs_genes <- class(slabel)
      cs_genes <- cs_genes[1]
      if(cs_genes == "data.frame" || cs_genes == "data.table" || cs_genes == "tbl_df"){
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        colnames(slabel)[genecolnum] <- "genelist"
        slabel <- slabel$genelist
      }
    }
    return(slabel)
  }
  
  if(class(result) == "list"){
    def_volcano <- function(x){
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if (length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot", i, "\n"))
            title1 = names(result)[i]
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = 2.5,
                                   labSize = 3.0,
                                   subtitle = NULL,
                                   xlim = c(10,-10),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot", i, "\n"))
          title1 = names(result)[i]
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = 2.5,
                                 labSize = 3.0,
                                 subtitle = NULL,
                                 xlim = c(10,-10),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
        }
      }
    }
    def_save_volcano <- function(x){
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if (length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
        }
        if (length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot", i, "\n"))
            title1 = names(result)[i]
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = 2.5,
                                   labSize = 3.0,
                                   subtitle = NULL,
                                   xlim = c(10,-10),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
            if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot", i, "\n"))
          title1 = names(result)[i]
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = 2.5,
                                 labSize = 3.0,
                                 subtitle = NULL,
                                 xlim = c(10,-10),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
        }
      }
    }
    
    def_title_volcano <- function(x){
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if (length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
            title1 = readline(prompt = "Please provide the title for plot", i, ": ")
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = 2.5,
                                   labSize = 3.0,
                                   subtitle = NULL,
                                   xlim = c(10,-10),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
          title1 = readline(prompt = "Please provide the title for plot", i, ": ")
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = 2.5,
                                 labSize = 3.0,
                                 subtitle = NULL,
                                 xlim = c(10,-10),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
        }
      }
    }
    def_title_save_volcano <- function(x){
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if(length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = 2.5,
                                     labSize = 3.0,
                                     subtitle = NULL,
                                     xlim = c(10,-10),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
            title1 = readline(prompt = "Please provide the title for plot", i, ": ")
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = 2.5,
                                   labSize = 3.0,
                                   subtitle = NULL,
                                   xlim = c(10,-10),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
            if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
          title1 = readline(prompt = "Please provide the title for plot", i, ": ")
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = 2.5,
                                 labSize = 3.0,
                                 subtitle = NULL,
                                 xlim = c(10,-10),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
        }
      }
    }
    
    ndef_volcano <- function(x){
      ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
      lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if(length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
        }
        if(length(x) == 1){
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot", i, "\n"))
            title1 = names(result)[i]
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = ps,
                                   labSize = ls,
                                   subtitle = NULL,
                                   xlim = c(lim,-abs(lim)),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot", i, "\n"))
          title1 = names(result)[i]
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = ps,
                                 labSize = ls,
                                 subtitle = NULL,
                                 xlim = c(lim,-abs(lim)),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
        }
      }
    }
    ndef_save_volcano <- function(x){
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      
      ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
      lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if(length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot", i, "\n"))
              title1 = names(result)[i]
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot", i, "\n"))
            title1 = names(result)[i]
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = ps,
                                   labSize = ls,
                                   subtitle = NULL,
                                   xlim = c(lim,-abs(lim)),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
            if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot", i, "\n"))
          title1 = names(result)[i]
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = ps,
                                 labSize = ls,
                                 subtitle = NULL,
                                 xlim = c(lim,-abs(lim)),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
        }
      }
    }
    
    ndef_title_volcano <- function(x){
      ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
      lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if(length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
            title1 = readline(prompt = "Please provide the title for plot", i, ": ")
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = ps,
                                   labSize = ls,
                                   subtitle = NULL,
                                   xlim = c(lim,-abs(lim)),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
          title1 = readline(prompt = "Please provide the title for plot", i, ": ")
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = ps,
                                 labSize = ls,
                                 subtitle = NULL,
                                 xlim = c(lim,-abs(lim)),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
        }
      }
    }
    ndef_title_save_volcano <- function(x){
      ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
      lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(x)){
        np = x[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
      
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      
      if (tolower(slab) == "yes" || tolower(slab) == "y"){
        if(length(x) >= 2){
          slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")
          if (tolower(slab2) == "yes" || tolower(slab2) == "y"){
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              slabel <- slab_func()
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
          if (tolower(slab2) == "no" || tolower(slab2) == "n"){
            slabel <- slab_func()
            for (l in 1:length(x)){
              i = x[l]
              cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
              title1 = readline(prompt = "Please provide the title for plot", i, ": ")
              res <- as.data.frame(result[[i]][[2]])
              plot <-EnhancedVolcano(res, 
                                     lab = res$hgnc_symbol,
                                     selectLab = slabel,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = title1,
                                     pCutoff = p,
                                     FCcutoff = lfc,
                                     pointSize = ps,
                                     labSize = ls,
                                     subtitle = NULL,
                                     xlim = c(lim,-abs(lim)),
                                     legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
              )
              print(plot)
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            }
          }
        }
        if(length(x) == 1){
          slabel <- slab_func()
          for (l in 1:length(x)){
            i = x[l]
            cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
            title1 = readline(prompt = "Please provide the title for plot", i, ": ")
            res <- as.data.frame(result[[i]][[2]])
            plot <-EnhancedVolcano(res, 
                                   lab = res$hgnc_symbol,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title1,
                                   pCutoff = p,
                                   FCcutoff = lfc,
                                   pointSize = ps,
                                   labSize = ls,
                                   subtitle = NULL,
                                   xlim = c(lim,-abs(lim)),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
            )
            print(plot)
            if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
            if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          }
        }
      }
      if(tolower(slab) == "no" || tolower(slab) == "no"){
        for (l in 1:length(x)){
          i = x[l]
          cat(paste("Now printing plot ", i, ". ", names(reesult)[i], "\n", sep = ""))
          title1 = readline(prompt = "Please provide the title for plot", i, ": ")
          res <- as.data.frame(result[[i]][[2]])
          plot <-EnhancedVolcano(res, 
                                 lab = res$hgnc_symbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = title1,
                                 pCutoff = p,
                                 FCcutoff = lfc,
                                 pointSize = ps,
                                 labSize = ls,
                                 subtitle = NULL,
                                 xlim = c(lim,-abs(lim)),
                                 legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
          )
          print(plot)
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste0(dir_n, "/", names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(paste0(names(result)[i], ".png", sep = ""), width = 8, height = 8, dpi = 500)}
        }
      }
    }
    
    p = as.numeric(readline(prompt = "Enter the desired adjusted p-value: "))
    lfc = as.numeric(readline(prompt = "Enter the desired LFC cutoff: "))
    
    all_res <- readline(prompt = "Do you want to plot all your results? (yes/no): ")
    if(tolower(all_res) == "yes" || tolower(all_res) == "y"){
      nplot <- as.numeric(seq(1, length(result)))
      def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
      p_title <- readline(prompt = "Do you want to provide your own plot titles? (yes/no): ")
      save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
      
      if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){def_title_save_volcano(nplot)}
      if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "no"){def_title_volcano(nplot)}
      if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "yes"){def_save_volcano(nplot)}
      if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "no"){def_volcano(nplot)}
      if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_volcano(nplot)}
      if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "no"){ndef_title_volcano(nplot)}
      if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "yes"){ndef_save_volcano(nplot)}
      if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "no"){ndef_volcano(nplot)}
    }
    
    if(tolower(all_res) == "no" || tolower(all_res) == "n"){
      message("The following are your provided results: ")
      for (np in 1:length(result)){cat(np, ". ", names(result)[np], "\n", sep = "")}
      num_plot <- as.numeric(readline(prompt = "Please state the number of plots required: "))
      
      if(num_plot >= 2){
        message("Please specify the result number for the plots in a comma separated format WITH ONE SPACE: e.g. 1, 2, 4")
        splot <- readline(prompt = "Please enter the result number of the plots required: ")
        splotn <- as.numeric(strsplit(splot, ", ")[[1]]) 
        def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
        p_title <- readline(prompt = "Do you want to provide your own plot titles? (yes/no): ")
        save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
        
        if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){def_title_save_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "no"){def_title_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "yes"){def_save_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "no"){def_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "no"){ndef_title_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "yes"){ndef_save_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "no"){ndef_volcano(splotn)}
      }
      
      if(num_plot == 1){
        splotn <- as.numeric(readline(prompt = "Please specify the result number for the plot required: "))
        def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
        p_title <- readline(prompt = "Do you want to provide your own plot titles? (yes/no): ")
        save_res <- readline(prompt = "Do you want to save your plots? (yes/no): ")
        
        if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){def_title_save_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "yes" && tolower(save_res) == "no"){def_title_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "yes"){def_save_volcano(splotn)}
        if(tolower(def_par) == "yes" && tolower(p_title) == "no" && tolower(save_res) == "no"){def_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "yes"){ndef_title_save_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "yes" && tolower(save_res) == "no"){ndef_title_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "yes"){ndef_save_volcano(splotn)}
        if(tolower(def_par) == "no" && tolower(p_title) == "no" && tolower(save_res) == "no"){ndef_volcano(splotn)}
      }
    }
  }
  if(class(result) == "DESeqResults"){
    dsr_def_volcano <- function(x){
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      if(tolower(slab) == "yes" || tolower(slab) == "y"){
        slabel <- slab_func()
        title1 = readline(prompt = "Please provide a title for the plot: ")
        res <- as.data.frame(x)
        plot <-EnhancedVolcano(res, 
                               lab = rownames(res),
                               selectLab = slabel,
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = title1,
                               pCutoff = p,
                               FCcutoff = lfc,
                               pointSize = 2.5,
                               labSize = 3.0,
                               subtitle = NULL,
                               xlim = c(10,-10),
                               legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
        )
        print(plot)
        save_res <- readline(prompt = "Do you want to save your plot? (yes/no): ")
        if(tolower(save_res) == "yes" || tolower(save_res) == "y"){ggsave(paste0(title1, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
      }
      if(tolower(slab) == "no" || tolower(slab) == "n"){
        title1 = readline(prompt = "Please provide a title for the plot: ")
        res <- as.data.frame(x)
        plot <-EnhancedVolcano(res, 
                               lab = rownames(res),
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = title1,
                               pCutoff = p,
                               FCcutoff = lfc,
                               pointSize = 2.5,
                               labSize = 3.0,
                               subtitle = NULL,
                               xlim = c(10,-10),
                               legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
        )
        print(plot)
        save_res <- readline(prompt = "Do you want to save your plot? (yes/no): ")
        if(tolower(save_res) == "yes" || tolower(save_res) == "y"){ggsave(paste0(title1, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
      }
    }
    dsr_ndef_volcano <- function(x){
      slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
      if(tolower(slab) == "yes" || tolower(slab) == "y"){
        slabel <- slab_func()
        title1 = readline(prompt = "Please provide a title for the plot: ")
        ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
        ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
        lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
        res <- as.data.frame(x)
        plot <-EnhancedVolcano(res, 
                               lab = rownames(res),
                               selectLab = slabel,
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = title1,
                               pCutoff = p,
                               FCcutoff = lfc,
                               pointSize = ps,
                               labSize = ls,
                               subtitle = NULL,
                               xlim = c(lim,-abs(lim)),
                               legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
        )
        print(plot)
        save_res <- readline(prompt = "Do you want to save your plot? (yes/no): ")
        if(tolower(save_res) == "yes" || tolower(save_res) == "y"){ggsave(paste0(title1, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
      }
      if(tolower(slab) == "no" || tolower(slab) == "n"){
        title1 = readline(prompt = "Please provide a title for the plot: ")
        ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
        ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.0: "))
        lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
        res <- as.data.frame(x)
        plot <-EnhancedVolcano(res, 
                               lab = rownames(res),
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = title1,
                               pCutoff = p,
                               FCcutoff = lfc,
                               pointSize = ps,
                               labSize = ls,
                               subtitle = NULL,
                               xlim = c(lim,-abs(lim)),
                               legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
        )
        print(plot)
        save_res <- readline(prompt = "Do you want to save your plot? (yes/no): ")
        if(tolower(save_res) == "yes" || tolower(save_res) == "y"){ggsave(paste0(title1, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
      }
    }
    
    def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
    
    if(tolower(def_par) == "yes" || tolower(def_par) == "y"){dsr_def_volcano(result)}
    if(tolower(def_par) == "no" || tolower(def_par) == "n"){dsr_ndef_volcano(result)}
  }
}

countplot <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }
  
  slab_func <- function(){
    message("Please ensure that the genes provided are in official gene symbol")
    message("Do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
    input = readline(prompt = "Enter D or V: ")
    if(tolower(input) == "d"){
      message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
      genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
      slabel <- strsplit(genes, ", ")[[1]]
    }
    else if(tolower(input) == 'v'){
      message("Please provide a dataframe or a vector containing the gene names :)")
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        slabel <- get(genes)
        
      }, error = function(e){
        message("Error: ", e)
        message("Please provide a valid dataframe or a vector of gene names: ")
        return(NULL)
      })
      
      cs_genes <- class(slabel)
      cs_genes <- cs_genes[1]
      if(cs_genes == "data.frame" || cs_genes == "data.table" || cs_genes == "tbl_df"){
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        colnames(slabel)[genecolnum] <- "geneid"
        slabel <- slabel$geneid
      }
    }
    return(slabel)
  }
  
  hugo <- readline(prompt = "Are the gene names provided in the expression data official HUGO gene symbols? (yes/no): ")
  plot_all <- readline(prompt = "Do you want to plot all the genes in the provided dataset? (yes/no): ")
  if(tolower(plot_all) == "yes" || tolower(plot_all) == "y"){genelist = rownames(assay(dds))}
  if(tolower(plot_all) == "no" || tolower(plot_all) == "n"){genelist <- slab_func()}
  
  if(tolower(hugo) == "yes" || tolower(hugo) == "y"){
    plot_n <- readline(prompt = "Do you want more than 1 graph/plot per picture? (yes/no): ")
    if(tolower(plot_n) == "yes" || tolower(plot_n) == "y"){
      cat("Please note that detailed parameters cannot be altered and you'll have to save it individually \n")
      proc <- readline(prompt = "Do you wish to change your mind? (yes/no): ")
      
      if(tolower(proc) == "no" || tolower(proc) == "n"){
        cat("Please state the number of rows and the number of plots per row\n")
        cat("For example, 2 rows and 3 plots per row will produce a 2x3 picture\n")
        row_n <- as.numeric(readline(prompt = "Please state the number of rows: "))
        column_n <- as.numeric(readline(prompt = "Please state the number of plot per row: "))
        par(mfrow = c(row_n, column_n))
        for (i in unique(1:length(genelist))){
          gp <- genelist[i]
          tryCatch({
            plotCounts(dds, gene = gp, intgroup = "condition")
          }, error = function(e) {
            cat(gp, "not found \n")
          })
        }
      }
    }
    
    def_par <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")
    save_res <- readline(prompt = "Do you want to save your results? (yes/no): ")
    if(tolower(def_par) == "yes" || tolower(def_par) == "y" && tolower(save_res) == "yes" || tolower(save_res) == "y"){
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      for (i in unique(1:length(genelist)))
      {
        gp <- genelist[i]
        tryCatch({
          gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
          cat("Now printing the count plot for gene", gp, "\n")
          ggraph = 
            ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
            geom_boxplot(width = 0.3) +
            geom_beeswarm(cex = 1.0, size = 1.0) +
            scale_y_log10(labels = scales::comma) +
            xlab("Condition")+
            ylab("Normalised Counts")+
            scale_color_brewer(palette = "Set1") +
            ggtitle(gp) +
            theme_pubclean() + 
            theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
                  axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                  axis.title = element_text(size = 16), axis.text = element_text(size = 15))
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(file = paste(dir_n, "/", gp, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(file = paste(gp, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          print(ggraph)
        }, error = function(e) {
          cat(gp, "not found \n")
        })
      }
    }
    if(tolower(def_par) == "no" || tolower(def_par) == "n" && tolower(save_res) == "yes" || tolower(save_res) == "y"){
      dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exist in your current directory")
        dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
        dir.create(dir_n)
      }
      box_w <- as.numeric(readline(prompt = "Please enter the width of the boxplot. Default = 0.5: "))
      labx = readline(prompt = "Please enter the x-axis label: ")
      laby = readline(prompt = "Please enter the y-axis label: ")
      t_size = as.numeric(readline(prompt = "Please enter the size of the title. Default = 18: "))
      a_size = as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
      at_size = as.numeric(readline(prompt = "Please enter the size of the axis tick labels. Default = 15: "))
      
      for (i in unique(1:length(genelist)))
      {
        gp <- genelist[i]
        tryCatch({
          gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
          cat("Now printing the count plot for gene", gp, "\n")
          ggraph = 
            ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
            geom_boxplot(width = box_w) +
            geom_beeswarm(cex = 1.0, size = 1.0) +
            scale_y_log10(labels = scales::comma) +
            xlab(labx)+
            ylab(laby)+
            scale_color_brewer(palette = "Set1") +
            ggtitle(gp) +
            theme_pubclean() + 
            theme(plot.title = element_text(hjust = 0.5, size = t_size), legend.position = 'none',
                  axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                  axis.title = element_text(size = a_size), axis.text = element_text(size = at_size))
          if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(file = paste(dir_n, "/", gp, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(file = paste(gp, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
          print(ggraph)
        }, error = function(e) {
          cat(gp, "not found \n")
        })
      }
    }
    
    if(tolower(def_par) == "yes" || tolower(def_par) == "y" && tolower(save_res) == "no" || tolower(save_res) == "n"){
      for (i in unique(1:length(genelist)))
      {
        gp <- genelist[i]
        tryCatch({
          gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
          cat("Now printing the count plot for gene", gp, "\n")
          ggraph = 
            ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
            geom_boxplot(width = 0.3) +
            geom_beeswarm(cex = 1.0, size = 1.0) +
            scale_y_log10(labels = scales::comma) +
            xlab("Condition")+
            ylab("Normalised Counts")+
            scale_color_brewer(palette = "Set1") +
            ggtitle(gp) +
            theme_pubclean() + 
            theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
                  axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                  axis.title = element_text(size = 16), axis.text = element_text(size = 15))
          print(ggraph)
        }, error = function(e) {
          cat(gp, "not found \n")
        })
      }
    }
    if(tolower(def_par) == "no" || tolower(def_par) == "n" && tolower(save_res) == "no" || tolower(save_res) == "n"){
      box_w <- as.numeric(readline(prompt = "Please enter the width of the boxplot. Default = 0.5: "))
      labx = readline(prompt = "Please enter the x-axis label: ")
      laby = readline(prompt = "Please enter the y-axis label: ")
      t_size = as.numeric(readline(prompt = "Please enter the size of the title. Default = 18: "))
      a_size = as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
      at_size = as.numeric(readline(prompt = "Please enter the size of the axis tick labels. Default = 15: "))
      
      for (i in unique(1:length(genelist)))
      {
        gp <- genelist[i]
        tryCatch({
          gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
          cat("Now printing the count plot for gene", gp, "\n")
          ggraph = 
            ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
            geom_boxplot(width = box_w) +
            geom_beeswarm(cex = 1.0, size = 1.0) +
            scale_y_log10(labels = scales::comma) +
            xlab(labx)+
            ylab(laby)+
            scale_color_brewer(palette = "Set1") +
            ggtitle(gp) +
            theme_pubclean() + 
            theme(plot.title = element_text(hjust = 0.5, size = t_size), legend.position = 'none',
                  axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                  axis.title = element_text(size = a_size), axis.text = element_text(size = at_size))
          print(ggraph)
        }, error = function(e) {
          cat(gp, "not found \n")
        })
      }
    }
  }
  
  if(tolower(hugo) == "no" || tolower(hugo) == "n"){
   resn <- readline(prompt = "Please enter the name of your results dataframe or list from the getresults function: ")
   tryCatch({
     res <- get(resn)
     
   }, error = function(e){
     message("Error: ", e)
     message("Please provide a valid dataframe or a list of results: ")
     return(NULL)
   })
   res <- res[[1]][[2]]
   
   if(tolower(plot_all) == "yes" || tolower(plot_all) == "y"){
     if(any(grepl("ENS", genelist))){genelist2 <- gsub("\\.\\d+$", "", genelist)}
     gene_index <- match(genelist2, res$geneid)
     genelist_hgnc <- res$hgnc_symbol[gene_index]
   }
   
   if(tolower(plot_all) == "no" || tolower(plot_all) == "n"){
     genelist_hgnc <- genelist
     gene_index <- match(genelist, res$hgnc_symbol)
     genelist <- res$geneid[gene_index]
     if(any(grepl("ENS", rownames(dds)))){rownames(dds) <- gsub("\\.\\d+$", "", rownames(dds))}
   }

   def_par <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")
   save_res <- readline(prompt = "Do you want to save your results? (yes/no): ")
   if(tolower(def_par) == "yes" || tolower(def_par) == "y" && tolower(save_res) == "yes" || tolower(save_res) == "y"){
     dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
     if(tolower(dir) == "yes" || tolower(dir) == "y"){
       message("Please make sure the name of your created folder does not exist in your current directory")
       dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
       dir.create(dir_n)
     }
     for (i in unique(1:length(genelist)))
     {
       gp <- genelist[i]
       gt <- genelist_hgnc[i]
       tryCatch({
         gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
         cat("Now printing the count plot for gene", gt, "\n")
         ggraph = 
           ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
           geom_boxplot(width = 0.3) +
           geom_beeswarm(cex = 1.0, size = 1.0) +
           scale_y_log10(labels = scales::comma) +
           xlab("Condition")+
           ylab("Normalised Counts")+
           scale_color_brewer(palette = "Set1") +
           ggtitle(gt) +
           theme_pubclean() + 
           theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
                 axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                 axis.title = element_text(size = 16), axis.text = element_text(size = 15))
         if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(file = paste(dir_n, "/", gt, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
         if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(file = paste(gt, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
         print(ggraph)
       }, error = function(e) {
         cat(gp, "not found \n")
       })
     }
   }
   if(tolower(def_par) == "no" || tolower(def_par) == "n" && tolower(save_res) == "yes" || tolower(save_res) == "y"){
     dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
     if(tolower(dir) == "yes" || tolower(dir) == "y"){
       message("Please make sure the name of your created folder does not exist in your current directory")
       dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
       dir.create(dir_n)
     }
     box_w <- as.numeric(readline(prompt = "Please enter the width of the boxplot. Default = 0.5: "))
     labx = readline(prompt = "Please enter the x-axis label: ")
     laby = readline(prompt = "Please enter the y-axis label: ")
     t_size = as.numeric(readline(prompt = "Please enter the size of the title. Default = 18: "))
     a_size = as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
     at_size = as.numeric(readline(prompt = "Please enter the size of the axis tick labels. Default = 15: "))
     
     for (i in unique(1:length(genelist)))
     {
       gp <- genelist[i]
       gt <- genelist_hgnc[i]
       tryCatch({
         gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
         cat("Now printing the count plot for gene", gt, "\n")
         ggraph = 
           ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
           geom_boxplot(width = box_w) +
           geom_beeswarm(cex = 1.0, size = 1.0) +
           scale_y_log10(labels = scales::comma) +
           xlab(labx)+
           ylab(laby)+
           scale_color_brewer(palette = "Set1") +
           ggtitle(gt) +
           theme_pubclean() + 
           theme(plot.title = element_text(hjust = 0.5, size = t_size), legend.position = 'none',
                 axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                 axis.title = element_text(size = a_size), axis.text = element_text(size = at_size))
         if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(file = paste(dir_n, "/", gt, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
         if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(file = paste(gt, ".png", sep = ""), width = 8, height = 8, dpi = 500)}
         print(ggraph)
       }, error = function(e) {
         cat(gp, "not found \n")
       })
     }
   }
   
   if(tolower(def_par) == "yes" || tolower(def_par) == "y" && tolower(save_res) == "no" || tolower(save_res) == "n"){
     for (i in unique(1:length(genelist)))
     {
       gp <- genelist[i]
       gt <- genelist_hgnc[i]
       tryCatch({
         gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
         cat("Now printing the count plot for gene", gt, "\n")
         ggraph = 
           ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
           geom_boxplot(width = 0.3) +
           geom_beeswarm(cex = 1.0, size = 1.0) +
           scale_y_log10(labels = scales::comma) +
           xlab("Condition")+
           ylab("Normalised Counts")+
           scale_color_brewer(palette = "Set1") +
           ggtitle(gt) +
           theme_pubclean() + 
           theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
                 axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                 axis.title = element_text(size = 16), axis.text = element_text(size = 15))
         print(ggraph)
       }, error = function(e) {
         cat(gp, "not found \n")
       })
     }
   }
   if(tolower(def_par) == "no" || tolower(def_par) == "n" && tolower(save_res) == "no" || tolower(save_res) == "n"){
     box_w <- as.numeric(readline(prompt = "Please enter the width of the boxplot. Default = 0.5: "))
     labx = readline(prompt = "Please enter the x-axis label: ")
     laby = readline(prompt = "Please enter the y-axis label: ")
     t_size = as.numeric(readline(prompt = "Please enter the size of the title. Default = 18: "))
     a_size = as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
     at_size = as.numeric(readline(prompt = "Please enter the size of the axis tick labels. Default = 15: "))
     
     for (i in unique(1:length(genelist)))
     {
       gp <- genelist[i]
       gt <- genelist_hgnc[i]
       tryCatch({
         gplot <- plotCounts(dds, gene= gp , intgroup = "condition", returnData = TRUE)
         cat("Now printing the count plot for gene", gt, "\n")
         ggraph = 
           ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
           geom_boxplot(width = box_w) +
           geom_beeswarm(cex = 1.0, size = 1.0) +
           scale_y_log10(labels = scales::comma) +
           xlab(labx)+
           ylab(laby)+
           scale_color_brewer(palette = "Set1") +
           ggtitle(gt) +
           theme_pubclean() + 
           theme(plot.title = element_text(hjust = 0.5, size = t_size), legend.position = 'none',
                 axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                 axis.title = element_text(size = a_size), axis.text = element_text(size = at_size))
         print(ggraph)
       }, error = function(e) {
         cat(gp, "not found \n")
       })
     }
   }
  }
}

datatransform <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }
  
  message("There are three types of transformation for downstream visualisations:")
  cat("1. Variance Stabilising Transformation (vst)\n")
  cat("2. Normalised Count Transformation (nct)\n")
  cat("3. Regularised Log Transformation(rlog)\n")
  message("Please check the DESeq2 vignette below for further details on the transformations: ")
  cat("http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html\n")
  
  tmtd <- readline(prompt = "Please state the desired method of transformation (vst/nct/rlog): ")
  
  if(tolower(tmtd) == "vst"){
    bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (Default = TRUE) (TRUE/FALSE): ")
    bt %<>% as.logical()
    tdds <- vst(dds, blind = bt)
  }
  
  if(tolower(tmtd) == "nct"){
    cat("Parameters for nct include function and pseudocount..\n")
    cat("Default function: log2\n")
    cat("Default pseudocount: 1\n")
    pcount <- as.numeric(readline(prompt = "Please state the desired pseudocount. Default = 1: "))
    tdds <- normTransform(dds, f = log2, pc = pcount)
  }
  
  if(tolower(tmtd) == "rlog"){
    bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (Default = TRUE) (TRUE/FALSE): ")
    bt %<>% as.logical()
    tdds <- rlog(dds, blind = bt)
  }
  
  pca_res = readline(prompt = "Do you want to plot the PCA Plot? (yes/no): ")
  if(tolower(pca_res) == "yes" || tolower(pca_res) == "y"){
    t <- readline(prompt = "Please provide the title of the plot: ")
    llab <- readline(prompt = "Please provide the legend title of the plot: ")
    plotPCA(tdds, intgroup = "condition")
    pcaplot <- plotPCA(tdds, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pcaplot, "percentVar"))
    graph <- ggplot(pcaplot, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed() +
      theme_minimal() +
      scale_fill_brewer() +
      labs(color = llab) +
      ggtitle(t) 
    print(graph)
  }
  
  hmap_samp = readline(prompt = "Do you want to plot the sample-to-sample distances heatmap? (yes/no): ")
  if(tolower(hmap_samp) == "yes" || tolower(hmap_samp) == "y"){
    sampleDists = dist(t(assay(tdds)))
    samplematrix = as.matrix(sampleDists)
    rownames(samplematrix) <- paste(tdds$condition, sep="-")
    colnames(samplematrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
    pheatmap(samplematrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors) %>% print()
  }
  save_res <- readline(prompt = "Do you want to save the transformed result? (yes/no): ")
  if(tolower(save_res) == "yes" || tolower(save_res) == "y"){
    hugo <- readline(prompt = "Are the genes in the provided expression dataset official HUGO gene symbols? (yes/no): ")
    if(tolower(hugo) == "yes" || tolower(hugo) == "y"){
      tres <- assay(tdds) %>% as.data.frame()
      tres$hgnc_symbol <- rownames(tdds)
      tres %<>% relocate(hgnc_symbol)
      rownames(tres) <- NULL
      filen <- readline(prompt = "Please name your transformed file: ")
      write_xlsx(tres, paste("./", filen, "_transformedcd.xlsx", sep = ""))
    }
    if(tolower(hugo) == "no" || tolower(hugo) == "n"){
      src <- readline(prompt = "Are the gene names provided in the expression data Ensembl Gene IDs? (yes/no): ")
      if(src == "yes" || src == "y"){
        srcType = "ensembl_gene_id"
      }
      if(src == "no" || src == "n"){
        message("Please refer to the specific attributte name in attrlist")
        srcType <- readline(prompt = "Please enter the specific attribute name: ")
      }

      if(srcType == "ensembl_gene_id"){
        tres <- assay(tdds) %>% as.data.frame()
        tres$geneid <- sapply(strsplit(rownames(tres),split="\\+"), "[",1)
        tres$geneid <- gsub("\\.\\d+$", "", tres$geneid)}
        v <- tres$geneid
        cat("Please wait.....\n")
        ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
      
        ## Make sure there was at least one mapping
        if( nrow(ID) < 1 ) top_n( "No IDs mapped successfully" )
        
        ## Drop empty duds
        k <- which( ID[,2] == "" )
        if( length(k) > 0 ) ID <- ID[-k,]
        colnames(ID)[1] <- "geneid"
        
        
        rem_res <- tres[!tres$geneid %in% ID[,1], ]
        rem_res %<>% dplyr::select(-geneid)
        rem_res$geneid <- rownames(rem_res)
        tres <- merge(tres, ID, by = "geneid", all.x = TRUE)
        tres %<>% dplyr::select(-geneid)
        tres %<>% relocate(hgnc_symbol)
        rownames(tres) <- NULL
        
        filen <- readline(prompt = "Please name your transformed file: ")
        write_xlsx(tres, paste("./", filen, "_transformedcd.xlsx", sep = ""))
        write_xlsx(rem_res, paste("./", filen, "_transformedcd_removedgenes.xlsx", sep = ""))
    }
  } 
  return(tdds)
}

cheatmap <- function(tdds){
  cs_dds <- class(tdds) == "DESeqTransform"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqTransform from the datatransform function")
    return(NULL)
  }
  
  slab_func <- function(){
    message("Do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
    input = readline(prompt = "Enter D or V: ")
    if(tolower(input) == "d"){
      message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
      genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
      slabel <- strsplit(genes, ", ")[[1]]
    }
    else if(tolower(input) == 'v'){
      message("Please provide a dataframe or a vector containing the gene names :)")
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        slabel <- get(genes)
        
      }, error = function(e){
        message("Error: ", e)
        message("Please provide a valid dataframe or a vector of gene names: ")
        return(NULL)
      })
      
      cs_genes <- class(slabel)
      cs_genes <- cs_genes[1]
      if(cs_genes == "data.frame" || cs_genes == "data.table" || cs_genes == "tbl_df"){
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        colnames(slabel)[genecolnum] <- "genelist"
        slabel <- slabel$genelist
      }
    }
    return(slabel)
  }
  
  genelist <- slab_func()
  hugo <- readline(prompt = "Are the genes in the provided expression data in official HUGO gene symbols? (yes/no): ")
  
  if(tolower(hugo) == "yes" || tolower(hugo) == "y"){
    llab <- readline(prompt = "Please provide the legend title: ")
    scn <- as.logical(readline(prompt = "Do you want to show the sample names? (TRUE/FALSE): "))
    ccol <- as.logical(readline(prompt = "Do you want to show the clustering between the samples? (TRUE/FALSE): "))
    
    sample = assay(tdds)[genelist,]
    sample <- (sample - rowMeans(sample))/rowSds(sample)
    collabel <- as.data.frame(colData(tdds)[c('condition')])
    colnames(collabel)[1] <- llab
    colour <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
    paletteLength <- 50
    myBreaks <- c(seq(min(sample), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(sample)/paletteLength, max(sample), length.out=floor(paletteLength/2)))
    pheatmap(sample, annotation_col = collabel, show_colnames = scn, cluster_cols = ccol, labels_row = genelist, color = colour, breaks = myBreaks) %>% print()
  }
  
  if(tolower(hugo) == "no" || tolower(hugo) == "n"){
    resn <- readline(prompt = "Please enter the name of your results dataframe or list from the getresults function: ")
    tryCatch({
      res <- get(resn)
      
    }, error = function(e){
      message("Error: ", e)
      message("Please provide a valid dataframe or a list of results: ")
      return(NULL)
    })
    res <- res[[1]][[2]]
    
    genelist_hgnc <- genelist
    gene_index <- match(genelist, res$hgnc_symbol)
    genelist <- res$geneid[gene_index]
    if(any(grepl("ENS", rownames(tdds)))){rownames(tdds) <- gsub("\\.\\d+$", "", rownames(tdds))}
    
    llab <- readline(prompt = "Please provide the legend title: ")
    scn <- as.logical(readline(prompt = "Do you want to show the sample names? (TRUE/FALSE): "))
    ccol <- as.logical(readline(prompt = "Do you want to show the clustering between the samples? (TRUE/FALSE): "))
    
    sample = assay(tdds)[genelist,]
    sample <- (sample - rowMeans(sample))/rowSds(sample)
    collabel <- as.data.frame(colData(tdds)[c('condition')])
    colnames(collabel)[1] <- llab
    colour <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
    paletteLength <- 50
    myBreaks <- c(seq(min(sample), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(sample)/paletteLength, max(sample), length.out=floor(paletteLength/2)))
    pheatmap(sample, annotation_col = collabel, show_colnames = scn, cluster_cols = ccol, labels_row = genelist, color = colour, breaks = myBreaks) %>% print()
  }
}