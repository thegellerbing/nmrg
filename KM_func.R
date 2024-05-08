genes2hgnc <- function(expr, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  expr %<>% as.data.frame()
  src <- readline(prompt = "Are the gene names provided in the expression data Ensembl Gene IDs? (yes/no): ")
  if(tolower(src) == "yes" || tolower(src) == "y"){
    expr[,1] <- gsub("\\.\\d+$", "", expr[[1]])
    srcType = "ensembl_gene_id"
  }
  if(tolower(src) == "no" || tolower(src) == "n"){
    attrlist = listAttributes(ensembl)
    message("Please proceed if you know the specific attribute in which to map your genes to")
    proc <- readline(prompt = "Do you want to proceed? (yes/no): ")
    if(tolower(proc) == "yes" || tolower(proc) == "y"){
      srcType = readline(prompt = "Please enter the specific attribute name in attrlist: ")
    }
    if(tolower(proc) == "no" || tolower(proc) == "n"){
      message("Please refer to the specific attribute in which to map your genes to")
      view(attrlist)
      return(attrlist)
    }
  }
  
  v <- expr[,1]
  
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top_n( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  colnames(ID)[1] <- colnames(expr)[1]
  expr <- merge(expr, ID, by = paste(colnames(expr)[1]))
  expr <- expr[, -1]
  expr %<>% relocate(hgnc_symbol)
  return(expr)
}

tidydata = function(){
  cat("Combining expression data with clinical data ....\n")
  expr_success <- FALSE
  while(!expr_success){
    n_expr <- readline(prompt = "Please enter the name of your expression data: ")
    tryCatch({
      expr <- get(n_expr)
      expr_success <- TRUE
    }, error = function(e){
      message("Error: Please provide a valid dataframe")
      expr_success <- FALSE
    })
  }
  
  md_success <- FALSE
  while(!md_success){
    n_md <- readline(prompt = "Please enter the name of your survival data: ")
    tryCatch({
      md <- get(n_md)
      md_success <- TRUE
    }, error = function(e){
      message("Error: Please provide a valid dataframe")
      md_success <- FALSE
    })
  }
  
  expr %<>% as.data.frame()
  md %<>% as.data.frame()
  
  message("For your query genes, do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
  input = readline(prompt = "Enter D or V (d/v): ")
  while(tolower(input) != "d" && tolower(input) != "v"){
    message("Answer not given, please try again :)")
    input = readline(prompt = "Enter D or V (d/v): ")
  }
  if(tolower(input) == "d"){
    message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
    genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
    genestostudy <- strsplit(genes, ", ")[[1]]
  }
  else if(tolower(input) == 'v'){
    message("Please provide a dataframe or a vector containing the gene names :)")
    gene_success <- FALSE
    while(!gene_success){
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        genestostudy <- get(genes)
        gene_success <- TRUE
      }, error = function(e){
        message("Error: Please provide a valid dataframe or a vector of gene names: ")
        gene_success <- FALSE
      })
    }
    
  cs_gts <- class(genestostudy)
  cs_gts <- cs_gts[1]
  if(cs_gts == "data.frame" || cs_gts == "data.table" || cs_gts == "tbl_df"){
    suppressWarnings({
      genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
      while(is.na(genecolnum) == TRUE || genecolnum <= 0 || genecolnum > ncol(genestostudy)){
        message("Error detected, please try again :)")
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        }
      colnames(genestostudy)[genecolnum] <- "geneid"
      genestostudy <- genestostudy$geneid
    })
    }
  }
  genestostudy <- unique(genestostudy)
  
  t_expr <- readline(prompt = "Does your expression data contain more than one column of gene identifiers? (yes/no): ")
  while(tolower(t_expr) != "yes" && tolower(t_expr) != "no"){
    message("Answer not given, please try again :)")
    t_expr <- readline(prompt = "Does your expression data contain more than one column of gene identifiers? (yes/no): ")
  }
  if(tolower(t_expr) == "yes" || tolower(t_expr) == "y"){
    message("If yes, columns that are not the official gene symbols have to be removed :)")
    suppressWarnings({
      hgnc_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      while(is.na(hgnc_ncol) == TRUE || hgnc_ncol <= 0 || hgnc_ncol > ncol(expr)){
        message("Error detected, please try again :)")
        hgnc_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
        }
    })
    colnames(expr)[hgnc_ncol] <- "hgnc"
    expr %<>% filter(expr$hgnc %in% genestostudy)
    genestostudy <- expr$hgnc
    rm_genes <- genestostudy[!grepl(paste0("\\b", paste(expr$hgnc, collapse = "\\b|\\b"), "\\b"), genestostudy)]
    if(length(rm_genes) > 0){
      for(i in 1:length(rm_genes)){cat(rm_genes[i], "was not found in your expression data \n")}
    }
    message("If there are more than one column of gene identifiers to remove, please enter it in a comma separated format WITH spaces")
    message("Example: 1, 3, 5")
    Sys.sleep(0.5)
    message("If there is only one, just key in the specific column number")
    message("Example: 1")
    suppressWarnings({
      gene_ncol <- readline(prompt = "Please state the column number in a comma separated format with spaces: ")
      gene_ncol <- as.numeric(strsplit(gene_ncol, ", ")[[1]])
      while(is.na(gene_ncol) == TRUE){
        message("Error detected, please try again :)")
        gene_ncol <- readline(prompt = "Please state the column number in a comma separated format with spaces: ")
        gene_ncol <- as.numeric(strsplit(gene_ncol, ", ")[[1]])
      }
    })
    expr <- expr[, -c(gene_ncol, hgnc_ncol)]
  }
  if(tolower(t_expr) == "no" || tolower(t_expr) == "n"){
    suppressWarnings({
      gene_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      while(is.na(gene_ncol) == TRUE || gene_ncol <= 0 || gene_ncol > ncol(expr)){
        message("Error detected, please try again :> ")
        gene_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
        }
    })
    colnames(expr)[gene_ncol] <- "hgnc"
    expr %<>% filter(expr$hgnc %in% genestostudy)
    rm_genes <- genestostudy[!grepl(paste0("\\b", paste(expr$hgnc, collapse = "\\b|\\b"), "\\b"), genestostudy)]
    if(length(rm_genes) > 0){
      for(i in 1:length(rm_genes)){cat(rm_genes[i], "was not found in your expression data \n")}
    }
    genestostudy <- expr$hgnc
    expr <- expr[, -gene_ncol]
  }

  expr <- as.data.frame(t(expr))
  Sys.sleep(0.5)
  
  view(md)
  suppressWarnings({
    idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
    while(is.na(idcol) == TRUE || idcol <= 0 || idcol > ncol(md)){
      message("Error detected, please try again :)")
      idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
      }
  })
  colnames(md)[idcol] <- 'id'
  
  suppressWarnings({
    surv_type <- readline(prompt = "Are you looking to test overall survival, progression free survival, or both? (os/pfs/both): ")
    while(is.na(surv_type) == TRUE || surv_type != "os" && surv_type != "pfs" && surv_type != "both"){
      message("Error detected, please re-enter your answer :)")
      surv_type <- readline(prompt = "Are you looking into OS, PFS, or both? (OS/PFS/both): ")
    }
  })
  
  if(surv_type == "os"){
    suppressWarnings({
      os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(md) || os_col == idcol){
        message("Error detected, please try again :)")
        os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      }
    })
    os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    while(tolower(os_val) != "yes" && tolower(os_val) != "no"){
      message("Answer not given, please try again :)")
      os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    }
    
    if(tolower(os_val) == "yes" || tolower(os_val) == "y"){
      cat("Converting overall survival to days and checking for NA values..\n")
      md %<>% as_tibble()
      colnames(md)[os_col] <- "os_months"
      md %<>% filter(md$os_months != "NA")
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      fd = function(x){x*30}
      days <- apply(md[,os_col], 2, fd)
      days %<>% as.data.frame()
      colnames(days) <- "os"
      md <- cbind(md, days)
      md %<>% as.data.frame()
    }
    
    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[os_col] <- "os"
      md %<>% filter(md$os != 'NA')
      md$os <- as.numeric(md$os)
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      md %<>% as.data.frame()
    }
    
    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "vitalstatus"
    vs <- readline(prompt = "Please provide the exact way the 'deceased' status is written in your provided survival dataframe: ")
    while(vs %in% md$vitalstatus == FALSE){
      message("Status provided not found, please try again :> ")
      vs <- readline(prompt = "Please provide the exact way the 'deceased' status is written in your provided survival dataframe: ")
    }
    md %<>%
      mutate(statusos = as.numeric(
        ifelse(
          vitalstatus == vs,
          '1',
          '0'
        )))
    
    md %<>% dplyr::select(id, os, vitalstatus, statusos)
    md <- cbind(md, expr2)
  } 
  
  if(surv_type == "pfs"){
    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(md) || pfs_col == idcol){
        message("Error detected, please try again :)")
        pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      }
    })
    pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    while(tolower(pfs_val) != "yes" && tolower(pfs_val) != "no"){
      message("Answer not given, please try again :)")
      pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    }
    
    if(tolower(pfs_val) == "yes" || tolower(pfs_val) == "y"){
      cat("Converting progression free survival to days and checking for NA values..\n")
      md %<>% as_tibble()
      colnames(md)[pfs_col] <- "pfs_months"
      md %<>% filter(md$pfs_months != "NA")
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      fd = function(x){x*30}
      days <- apply(md[,pfs_col], 2, fd)
      days %<>% as.data.frame()
      colnames(days) <- "pfs"
      md <- cbind(md, days)
      md %<>% as.data.frame()
    }
    
    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[pfs_col] <- "pfs"
      md %<>% filter(md$os != 'NA')
      md$pfs <- as.numeric(md$pfs)
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      md %<>% as.data.frame()
    }
    q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "pfsstatus"
    vs <- readline(prompt = "Please provide the exact way the 'progressed' status is written in your provided survival dataframe: ")
    while(vs %in% md$pfsstatus == FALSE){
      message("Status provided not found, please try again :> ")
      vs <- readline(prompt = "Please provide the exact way the 'progressed' status is written in your provided survival dataframe: ")
    }
    md %<>%
      mutate(statuspfs = as.numeric(
        ifelse(
          pfsstatus == vs,
          '1',
          '0'
        )))
    
    md %<>% dplyr::select(id, os, pfsstatus, statuspfs)
    md <- cbind(md, expr2)
  }
  
  if(surv_type == "both"){
    suppressWarnings({
      os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(md) || os_col == idcol){
        message("Error detected, please try again :)")
        os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      }
    })
    os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    while(tolower(os_val) != "yes" && tolower(os_val) != "no"){
      message("Answer not given, please try again :)")
      os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    }
    
    if(tolower(os_val) == "yes" || tolower(os_val) == "y"){
      cat("Converting overall survival to days and checking for NA values..\n")
      md %<>% as_tibble()
      colnames(md)[os_col] <- "os_months"
      md %<>% filter(md$os_months != "NA")
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      fd = function(x){x*30}
      days <- apply(md[,os_col], 2, fd)
      days %<>% as.data.frame()
      colnames(days) <- "os"
      md <- cbind(md, days)
      md %<>% as.data.frame()
    }
    
    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[os_col] <- "os"
      md %<>% filter(md$os != 'NA')
      md$os <- as.numeric(md$os)
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      md %<>% as.data.frame()
    }
    
    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "vitalstatus"
    vs <- readline(prompt = "Please provide the exact way the 'deceased' status is written in your provided survival dataframe: ")
    while(vs %in% md$vitalstatus == FALSE){
      message("Status provided not found, please try again :> ")
      vs <- readline(prompt = "Please provide the exact way the 'deceased' status is written in your provided survival dataframe: ")
    }
    md %<>%
      mutate(statusos = as.numeric(
        ifelse(
          vitalstatus == vs,
          '1',
          '0'
        )))
    
    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(md) || pfs_col == idcol){
        message("Error detected, please try again :)")
        pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      }
    })
    pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    while(tolower(pfs_val) != "yes" && tolower(pfs_val) != "no"){
      message("Answer not given, please try again :)")
      pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    }
    
    if(tolower(pfs_val) == "yes" || tolower(pfs_val) == "y"){
      cat("Converting progression free survival to days and checking for NA values..\n")
      md %<>% as_tibble()
      colnames(md)[pfs_col] <- "pfs_months"
      md %<>% filter(md$pfs_months != "NA")
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      fd = function(x){x*30}
      days2 <- apply(md[,pfs_col], 2, fd)
      days2 %<>% as.data.frame()
      colnames(days2) <- "pfs"
      md <- cbind(md, days2)
      md %<>% as.data.frame()
    }
    
    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[pfs_col] <- "pfs"
      md %<>% filter(md$os != 'NA')
      md$pfs <- as.numeric(md$pfs)
      expr %>% filter(rownames(expr) %in% md$id) %>% as.data.frame() -> expr2
      rm_expr <- setdiff(rownames(expr), rownames(expr2))
      md %<>% as.data.frame()
    }
    w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(w) == TRUE || w <= 0 || w > ncol(md) || w == idcol || w == os_col){
        message("Error detected, please try again :)")
        w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[w] <- "pfsstatus"
    vs <- readline(prompt = "Please provide the exact way the 'progressed' status is written in your provided survival dataframe: ")
    while(vs %in% md$pfsstatus == FALSE){
      message("Status provided not found, please try again :> ")
      vs <- readline(prompt = "Please provide the exact way the 'progressed' status is written in your provided survival dataframe: ")
    }
    md %<>%
      mutate(statuspfs = as.numeric(
        ifelse(
          pfsstatus == vs,
          '1',
          '0'
        )))
    md %<>% dplyr::select(id, os, vitalstatus, statusos, pfs, pfsstatus, statuspfs)
    md <- cbind(md, expr2)
  }
  
  if(!identical (rownames(expr2), md$id)) {
    md %<>% filter(md$id %in% rownames(expr2))
    expr2 %<>% filter(rownames(expr2) %in% md$id)
  }
  
  if(all(md$os != 'NA')) {cat("Data is clean \n")}
  
  tryCatch(
    stopifnot(all(rownames(expr2) %in% md$id)),
    error = function(e){
      message("Samples of the expression and metadata aren't equal")
    }
  )
  
  for (i in (1:length(colnames(expr2)))){
      x = colnames(expr2)[i]
      y = genestostudy[i]
      colnum <- match(x, colnames(md))
      names(md)[colnum] <- y
    }

  colnames(md) <- gsub("-", "", colnames(md))
  
  for(i in 1:length(rm_expr)){
    x = rm_expr[i]
    cat("Sample", x, "is not found in your clinical data \n")
  }
  
  cat("Congrats :> It is safe to proceed :> \n")
  rownames(md) <- NULL
  return(md)
}

uniKM = function(kmdf){
  cat("Executing Kaplan Meier Analysis .... \n")
  kmdf %<>% as.data.table()
  
  unicolumn = c('Gene', 'Coefficient', 'Hazard Ratio', '95% CI', 'P-value', 'PH Check')
  kmtable = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
  colnames(kmtable) = unicolumn
  
  message("For your query genes, do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
  input = readline(prompt = "Enter D or V (d/v): ")
  while(tolower(input) != "d" && tolower(input) != "v"){
    message("Answer not given, please try again :)")
    input = readline(prompt = "Enter D or V (d/v): ")
  }
  if(tolower(input) == "d"){
    message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
    genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
    genestostudy <- strsplit(genes, ", ")[[1]]
  }
  else if(tolower(input) == 'v'){
    message("Please provide a dataframe or a vector containing the gene names :)")
    gene_success <- FALSE
    while(!gene_success){
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        genestostudy <- get(genes)
        gene_success <- TRUE
      }, error = function(e){
        message("Error: Please provide a valid dataframe or a vector of gene names: ")
        gene_success <- FALSE
      })
    }
    
    cs_gts <- class(genestostudy)
    cs_gts <- cs_gts[1]
    if(cs_gts == "data.frame" || cs_gts == "data.table" || cs_gts == "tbl_df"){
      suppressWarnings({
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        while(is.na(genecolnum) == TRUE || genecolnum <= 0 || genecolnum > ncol(genestostudy)){
          message("Error detected, please try again :)")
          genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        }
        colnames(genestostudy)[genecolnum] <- "geneid"
        genestostudy <- genestostudy$geneid
      })
    }
  }
  genestostudy <- gsub("-", "", genestostudy)
  genestostudy <- genestostudy[genestostudy %in% colnames(kmdf)]
  
  default_theme <- function() {
    theme_survminer() %+replace%
      theme(plot.title=element_text(hjust=0.5, size = 18),legend.box = "vertical")+
      theme(plot.subtitle = element_text(hjust = 0.5, vjust = -75))}
  
  kmplot <- function(df){
    graph = ggsurvplot(fit1, data=df, 
                       break.time.by = 500, 
                       xlab = "Time (Days)", ylab = "Survival Probability", 
                       palette = c( "steelblue2", "firebrick1"),  
                       title = (paste(x, " Survival Curve", sep = "")), ggtheme = default_theme(),
                       legend.labs = c(paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = ""),
                                       paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = "")),
                       legend.title = "Median Survival", legend = c(0.65,0.75),
                       font.main = c(18,"bold"), font.legend = c(16),
                       font.x = c(16, "bold"), font.y = c(16, "bold"), font.tickslab = c(15)) + labs(subtitle = pvalue)
    print(graph)
    return(graph)
  }
  
  kmplot_pfs <- function(df){
    graph = ggsurvplot(fit1, data=df, 
                       break.time.by = 500, 
                       xlab = "Time (Days)", ylab = "Survival Probability", 
                       palette = c( "steelblue2", "firebrick1"),  
                       title = (paste(x, " Progression Free Survival", sep = "")), ggtheme = default_theme(),
                       legend.labs = c(paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = ""),
                                       paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = "")),
                       legend.title = "Median PFS", legend = c(0.65,0.75),
                       font.main = c(18,"bold"), font.legend = c(16),
                       font.x = c(16, "bold"), font.y = c(16, "bold"), font.tickslab = c(15)) + labs(subtitle = pvalue)
    print(graph)
    return(graph)
  }
  
  kmplot_par <- function(df){
    graph = ggsurvplot(fit1, data=df, 
                       break.time.by = 500, 
                       xlab = "Time (Days)", ylab = "Survival Probability", 
                       palette = c( "steelblue2", "firebrick1"),  
                       title = (paste(x, " Survival Curve", sep = "")), ggtheme = custom_theme(),
                       legend.labs = c(paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = ""),
                                       paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = "")), 
                       legend.title = "Median Survival", legend = lpos,
                       font.main = c(msize,"bold"), font.legend = c(lsize), font.subtitle = c(psize),
                       font.x = c(axsize, "bold"), font.y = c(axsize, "bold"), font.tickslab = c(15)) + labs(subtitle = pvalue)
    print(graph)
    return(graph)
  }
  
  kmplot_par_pfs <- function(df){
    graph = ggsurvplot(fit1, data=df, 
                       break.time.by = 500, 
                       xlab = "Time (Days)", ylab = "Survival Probability", 
                       palette = c( "steelblue2", "firebrick1"),  
                       title = (paste(x, " Progression Free Survival", sep = "")), ggtheme = custom_theme(),
                       legend.labs = c(paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = ""),
                                       paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = "")), 
                       legend.title = "Median PFS", legend = lpos,
                       font.main = c(msize,"bold"), font.legend = c(lsize), font.subtitle = c(psize),
                       font.x = c(axsize, "bold"), font.y = c(axsize, "bold"), font.tickslab = c(15)) + labs(subtitle = pvalue)
    print(graph)
    return(graph)
  }
  
  message("Do you want to set the cutoff based off the median, quartile, custom %, z-score, or individual optimal cutoffs?")
  cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
  while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal"){
    message("Answer not given, please try again :)")
    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
  }
  
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })
  
  plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  while(tolower(plot) != "yes" && tolower(plot) != "no"){
    message("Answer not given, please try again :>")
    plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  }
  if(tolower(plot) == "yes" || tolower(plot) == "y"){
    if(length(genestostudy) > 10){
      plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      while(tolower(plot_sig) != "yes" && tolower(plot_sig) != "no" && tolower(plot_sig) != "y" && tolower(plot_sig) != "n"){
        message("Answer not given, please try again :>")
        plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      }
    } else{
      plot_sig <- "no"
    }
    
    plot_type <- readline(prompt = "Are you plotting overall survival curves or PFS curves? (osc/pfs): ")
    while(tolower(plot_type) != "osc" && tolower(plot_type) != "pfs"){
      message("Answer not given, please try again :>")
      plot_type <- readline(prompt = "Are you plotting overall survival curves or PFS curves? (osc/pfs): ")
    }
    
    def_par <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")
    while(tolower(def_par) != "yes" && tolower(def_par) != "y" && tolower(def_par) != "no" && tolower(def_par) != "n"){
      message("Answer not given, please try again :>")
      def_par <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")
    }
    
    if(tolower(def_par) == "no" || tolower(def_par) == "no"){
      aspect <- readline(prompt = "Which aspect(s) of the graph(s) do you want to alter? (all/legends/title and axis/pvalue): ")
      suppressWarnings({
        while(is.na(aspect) == TRUE || tolower(aspect) != "all" && tolower(aspect) != "legends" && tolower(aspect) != "title and axis" && tolower(aspect) != "pvalue"){
          message("Error detected, please try again :)")
          aspect <- readline(prompt = "Which aspect of the graph do you want to alter? (all/legends/title and axis/pvalue): ")
        }
      })
      
      if(tolower(aspect) == "all"){
        cat("You are required to enter the x and y coordinates for the position of the legends \n")
        cat("Please enter the comma-separated values WITH SPACES in the following format:0.65, 0.75\n")
        suppressWarnings({
          lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
          lpos <- strsplit(lpos, ", ")[[1]]
          while(length(lpos) != 2){
            message("Coordinate provided is incompatible, please try again :)")
            lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
            lpos <- strsplit(lpos, ", ")[[1]]
          }
          lpos %<>% as.numeric()
          while(is.na(lpos) == TRUE || lpos[[1]] == 0 || lpos[[1]] > 1 || lpos[[2]] == 0 || lpos[[2]] > 1){
            message("Coordinate provided is incompatible, please try again :)")
            lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
          }
        })
        lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. Default = 16: "))
        suppressWarnings({
          while(is.na(lsize) == TRUE){
            message("Error detected, please try again :)")
            lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. Default = 16: "))
          }
        })
        msize <- as.numeric(readline(prompt = "Please enter the size of the main title. Default = 18: "))
        suppressWarnings({
          while(is.na(msize) == TRUE){
            message("Error detected, please try again :)")
            msize <- as.numeric(readline(prompt = "Please enter the size of the main title. Default = 18: "))
          }
        })
        
        axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
        suppressWarnings({
          while(is.na(axsize) == TRUE){
            message("Error detected, please try again :)")
            axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
          }
        })
        
        psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. Default = 16: "))
        suppressWarnings({
          while(is.na(psize) == TRUE){
            message("Error detected, please try again :)")
            psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. Default = 16: "))
          }
        })
        pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
        suppressWarnings({
          while(is.na(pxc) == TRUE){
            message("Error detected, please try again :)")
            pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
          }
        })
        cat("For the y-coordinate, the higher the value, the lower the position of the label.")
        pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. Default = 75: "))
        suppressWarnings({
          while(is.na(pxc) == TRUE){
            message("Error detected, please try again :)")
            pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. Default = 75: "))
          }
        })
      }
      if(tolower(aspect) == "legends"){
        cat("You are required to enter the x and y coordinates for the position of the legends \n")
        cat("Please enter the comma-separated values WITH SPACES in the following format:0.65, 0.75\n")
        suppressWarnings({
          lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
          lpos <- strsplit(lpos, ", ")[[1]]
          while(length(lpos) != 2){
            message("Coordinate provided is incompatible, please try again :)")
            lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
            lpos <- strsplit(lpos, ", ")[[1]]
          }
          lpos %<>% as.numeric()
          while(is.na(lpos) == TRUE || lpos[[1]] == 0 || lpos[[1]] > 1 || lpos[[2]] == 0 || lpos[[2]] > 1){
            message("Coordinate provided is incompatible, please try again :)")
            lpos <- readline(prompt="Please enter the x and y coordinates between 0 and 1 for the legend labels: ")
          }
        })
        lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. Default = 16: "))
        suppressWarnings({
          while(is.na(lsize) == TRUE){
            message("Error detected, please try again :)")
            lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. Default = 16: "))
          }
        })
        msize <- as.numeric(18)
        axsize <- as.numeric(16)
        psize <- as.numeric(16)
        pxc <- as.numeric(0.5)
        pyc <- as.numeric(75)
      }
      if(tolower(aspect) == "title and axis"){
        msize <- as.numeric(readline(prompt = "Please enter the size of the main title. Default = 18: "))
        suppressWarnings({
          while(is.na(msize) == TRUE){
            message("Error detected, please try again :)")
            msize <- as.numeric(readline(prompt = "Please enter the size of the main title. Default = 18: "))
          }
        })
        
        axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
        suppressWarnings({
          while(is.na(axsize) == TRUE){
            message("Error detected, please try again :)")
            axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
          }
        })
        
        psize <- as.numeric(16)
        pxc <- as.numeric(0.5)
        pyc <- as.numeric(75)
        
        lsize <- as.numeric(16)
        lpos <- "0.65, 0.75"
        lpos <- strsplit(lpos, ", ")[[1]]
        lpos %<>% as.numeric()
      }
      if(tolower(aspect) == "pvalue"){
        msize <- as.numeric(18)
        axsize <- as.numeric(16)
        lsize <- as.numeric(16)
        lpos <- "0.65, 0.75"
        lpos <- strsplit(lpos, ", ")[[1]]
        lpos %<>% as.numeric()
        
        psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. Default = 16: "))
        suppressWarnings({
          while(is.na(psize) == TRUE){
            message("Error detected, please try again :)")
            psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. Default = 16: "))
          }
        })
        pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
        suppressWarnings({
          while(is.na(pxc) == TRUE){
            message("Error detected, please try again :)")
            pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
          }
        })
        cat("For the y-coordinate, the higher the value, the lower the position of the label.")
        pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. Default = 75: "))
        suppressWarnings({
          while(is.na(pxc) == TRUE){
            message("Error detected, please try again :)")
            pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. Default = 75: "))
          }
        })
      }
      
      custom_theme <- function() {
        theme_survminer() %+replace%
          theme(plot.title=element_text(hjust=0.5),legend.box = "vertical")+
          theme(plot.subtitle = element_text(hjust = pxc, vjust = -pyc))}
    }
    
    if(length(genestostudy) > 10){
      dir <- readline(prompt = "Do you want to create a new folder for your plots? (yes/no): ")
      while(tolower(dir) != "yes" && tolower(dir) != "y" && tolower(dir) != "no" && tolower(dir) != "n"){
        message("Answer not given, please try again :>")
        dir <- readline(prompt = "Do you want to create a new folder for your plots & results? (yes/no): ")
      }
      
      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exists in your current directory")
        success <- FALSE
        while(!success){
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          if(dir.exists(dir_n)){
            message("Folder already exists. Please enter a different name :)")
          } else{
            dir.create(dir_n)
            success <- TRUE
          }
        }
      }
    } else{
      dir <- "no"
    }
  }
  
  plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  while(tolower(plotsr) != "yes" && tolower(plotsr) != "no"){
    message("Answer not given, please try again :>")
    plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  }
  if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
    if(!file.exists("PH Check")){
      dir.create("PH Check")
    }
    if(length(genestostudy) > 10){
      plotsr_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      while(tolower(plotsr_sig) != "yes" && tolower(plotsr_sig) != "no" && tolower(plotsr_sig) != "y" && tolower(plotsr_sig) != "n"){
        message("Answer not given, please try again :>")
        plotsr_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      }
    } else{
      plotsr_sig <- "no"
    }
  }
  
  if(tolower(cutoff) == "median"){
    if(plot_type == "osc"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(plot_type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}

    for (y in unique(1:length(genestostudy)))
    {
      z = genestostudy[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {kmdf %>% dplyr::select(all_of(z)) %>% unlist() -> c} 
      mdn = median(c)
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          c >= mdn,
          'High',
          'Low'
        )))}
    rm_df <- kmdf %>% dplyr::select(where(~all(. == "High")))
    rm_genes <- colnames(rm_df)
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    Sys.sleep(0.5)

    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
    }
    genestostudy <- gsub("Expression", "", gtslist)

    for (i in unique(1:length(genestostudy))) {
      x = genestostudy[i]
      y = gtslist[i]
      
      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          unicox = coxph(f2, data = kmdf)
          zph = cox.zph(unicox)
          phcheck = zph$table[1,3]
          ucoef = summary(unicox)$coef[1,1]
          uhaz = summary(unicox)$coef[1,2]
          upvalue = summary(unicox)$coef[1,5]
          ci = confint(unicox, level = 0.95)
          ci <- exp(ci[1,])
          ci <- round(ci, digits = 2)
          uci <- paste(ci[1], '~', ci[2], sep = " ")
          unidf2 = data.frame(x, ucoef, uhaz, uci, upvalue, phcheck)
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)
          
          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            slow = round(unname(summary(fit1)$table[,'median'][1]),0)
            shigh = round(unname(summary(fit1)$table[,'median'][2]),0)
            nlow = unname(summary(fit1)$table[,'records'][1])
            nhigh = unname(summary(fit1)$table[,'records'][2])
            pval <- surv_pvalue(fit1)
            pval <- pval[[4]]
            pval <- gsub("p = ", "", pval)
            pvalue <- bquote(italic(p) == .(pval))
            if(tolower(plot_sig) == "yes"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
                if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
                if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
                if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
                
                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }
            
            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
              if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
              if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
              if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
              
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              dev.off()
            }
          }
        } 
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }
  
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom"){
    if(tolower(cutoff) == "custom"){
      qval_low <- as.numeric(readline(prompt = "Please enter the low % cutoff between 0 and 1: "))
      qval_high <- as.numeric(readline(prompt = "Please enter the high % cutoff between 0 and 1: "))
    }
    
    for (y in unique(1:length(genestostudy)))
    {
      z = genestostudy[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(all_of(z)) %>% unlist() -> c
      } 
      if(tolower(cutoff) == "quartile"){quartile <- quantile(c, probs = c(0.25, 0.75))}
      if(tolower(cutoff) == "custom"){quartile <- quantile(c, probs = c(qval_low, qval_high))}
      
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          c >= quartile[[2]],
          'High',
          ifelse(
            c <= quartile[[1]], 
            "Low",
            "Medium"
          )
        )))
    }
    rm_df <- kmdf %>% dplyr::select(where(~all(. == "High")))
    rm_genes <- colnames(rm_df)
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    Sys.sleep(0.5)
    
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
    }
    genestostudy <- gsub("Expression", "", gtslist)
    
    for (i in unique(1:length(genestostudy))) {
      x = genestostudy[i]
      y = gtslist[i]
      
      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          kmdf %>% dplyr::select(id, os, vitalstatus, status, all_of(y)) %>% as.data.frame() -> kmdf2
          kmdf2 %<>% filter(kmdf2[,5] != "Medium")
          if(plot_type == "osc"){survdf <- Surv(time = kmdf2$os, event = kmdf2$statusos)}
          if(plot_type == "pfs"){survdf <- Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)}
          unicox = coxph(f2, data = kmdf2)
          zph = cox.zph(unicox)
          phcheck = zph$table[1,3]
          ucoef = summary(unicox)$coef[1,1]
          uhaz = summary(unicox)$coef[1,2]
          upvalue = summary(unicox)$coef[1,5]
          ci = confint(unicox, level = 0.95)
          ci <- exp(ci[1,])
          ci <- round(ci, digits = 2)
          uci <- paste(ci[1], '~', ci[2], sep = " ")
          unidf2 = data.frame(x, ucoef, uhaz, uci, upvalue, phcheck)
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)
          
          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            fit1 = surv_fit(f2, data = kmdf2)
            slow = round(unname(summary(fit1)$table[,'median'][1]),0)
            shigh = round(unname(summary(fit1)$table[,'median'][2]),0)
            nlow = unname(summary(fit1)$table[,'records'][1])
            nhigh = unname(summary(fit1)$table[,'records'][2])
            pval <- surv_pvalue(fit1)
            pval <- pval[[4]]
            pval <- gsub("p = ", "", pval)
            pvalue <- bquote(italic(p) == .(pval))
            if(tolower(plot_sig) == "yes"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
                if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
                if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
                if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
                
                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }
            
            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
              if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
              if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
              if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
              
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }            
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              dev.off()
            }
          }
        } 
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }
  
  if(tolower(cutoff) == "zscore"){
    zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))
    
    for (y in unique(1:length(genestostudy))){
      z = genestostudy[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(all_of(z)) %>% unlist() -> c
      } 
      zsc = scale(c)
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          zsc >= zcut,
          'High',
          ifelse(
            zsc <= -abs(zcut), 
            "Low",
            "Medium"
          )
        )))
    }
    
    rm_df <- kmdf %>% dplyr::select(where(~all(. == "High")))
    rm_genes <- colnames(rm_df)
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    Sys.sleep(0.5)
    
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist){
      fctcheck <- fct_count(kmdf[[o]])
      if(length(rownames(fctcheck)) <3){kmdf %<>% dplyr::select(-o)}
      if(length(rownames(fctcheck))== 3){kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')}
    }
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    genestostudy <- gsub("Expression", "", gtslist)
    
    for (i in unique(1:length(genestostudy))) {
      x = genestostudy[i]
      y = gtslist[i]
      
      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          kmdf %>% as.data.frame() -> kmdf2
          kmdf2 %<>% dplyr::select(id, os, vitalstatus, status, all_of(y))
          kmdf2 %<>% filter(kmdf2[,5] != "Medium")
          if(plot_type == "osc"){survdf <- Surv(time = kmdf2$os, event = kmdf2$statusos)}
          if(plot_type == "pfs"){survdf <- Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)}
          unicox = coxph(f2, data = kmdf2)
          zph = cox.zph(unicox)
          phcheck = zph$table[1,3]
          ucoef = summary(unicox)$coef[1,1]
          uhaz = summary(unicox)$coef[1,2]
          ucoef = summary(unicox)$coef[1,2]
          upvalue = summary(unicox)$coef[1,5]
          ci = confint(unicox, level = 0.95)
          ci <- exp(ci[1,])
          ci <- round(ci, digits = 2)
          uci <- paste(ci[1], '~', ci[2], sep = " ")
          unidf2 = data.frame(x, ucoef, uhaz, uci, upvalue, phcheck)
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)
          
          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            fit1 = surv_fit(f2, data = kmdf2)
            slow = round(unname(summary(fit1)$table[,'median'][1]),0)
            shigh = round(unname(summary(fit1)$table[,'median'][2]),0)
            nlow = unname(summary(fit1)$table[,'records'][1])
            nhigh = unname(summary(fit1)$table[,'records'][2])
            pval <- surv_pvalue(fit1)
            pval <- pval[[4]]
            pval <- gsub("p = ", "", pval)
            pvalue <- bquote(italic(p) == .(pval))
            if(tolower(plot_sig) == "yes"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
                if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
                if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
                if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
                
                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }
            
            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
              if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
              if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
              if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
              
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }            
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              dev.off()
            }
          }
        } 
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }
  
  if(tolower(cutoff) == "optimal"){
    if(plot_type == "osc"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(plot_type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    cutoffcol = c("Gene", "Cutoff")
    cuttable = data.frame(matrix(nrow=0, ncol = length(cutoffcol)))
    colnames(cuttable) = cutoffcol
    
    for (y in unique(1:length(genestostudy))){
      z = genestostudy[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(all_of(z)) %>% unlist() -> c
        fc <- as.formula(paste('survdf ~', paste(z)))
        tryCatch({
          mscutoff <- maxstat.test(fc, data = kmdf, smethod = "LogRank")
          o_cutoff <- mscutoff$estimate %>% as.double()
          cutdf <- data.frame(z, o_cutoff)
          cuttable <- rbind(cuttable, cutdf)
          
          kmdf %<>% mutate("{name}" := as.factor(
            ifelse(
              c >= o_cutoff,
              'High',
              'Low'
            )))
        }, error = function(e) {
          cat(z, "contains non-zero cell counts within the specified proportion range, removed from the analysis \n")
        })
      }
    }
    rm_df <- kmdf %>% dplyr::select(where(~all(. == "High")))
    rm_genes <- colnames(rm_df)
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    Sys.sleep(0.5)
   
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
    }
    genestostudy <- gsub("Expression", "", gtslist)
    
    for (i in unique(1:length(genestostudy))) {
      x = genestostudy[i]
      y = gtslist[i]
      
      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          unicox = coxph(f2, data = kmdf)
          zph = cox.zph(unicox)
          phcheck = zph$table[1,3]
          ucoef = summary(unicox)$coef[1,1]
          uhaz = summary(unicox)$coef[1,2]
          upvalue = summary(unicox)$coef[1,5]
          ci = confint(unicox, level = 0.95)
          ci <- exp(ci[1,])
          ci <- round(ci, digits = 2)
          uci <- paste(ci[1], '~', ci[2], sep = " ")
          unidf2 = data.frame(x, ucoef, uhaz, uci, upvalue, phcheck)
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)
          
          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            slow = round(unname(summary(fit1)$table[,'median'][1]),0)
            shigh = round(unname(summary(fit1)$table[,'median'][2]),0)
            nlow = unname(summary(fit1)$table[,'records'][1])
            nhigh = unname(summary(fit1)$table[,'records'][2])
            pval <- surv_pvalue(fit1)
            pval <- pval[[4]]
            pval <- gsub("p = ", "", pval)
            pvalue <- bquote(italic(p) == .(pval))
            if(tolower(plot_sig) == "yes"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
                if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
                if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
                if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
                
                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }
            
            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              if(tolower(def_par) == "yes" && tolower(plot_type) == "osc"){graph <- kmplot(kmdf)}
              if(tolower(def_par) == "no"  && tolower(plot_type) == "osc"){graph <- kmplot_par(kmdf)}
              if(tolower(def_par) == "yes" && tolower(plot_type) == "pfs"){graph <- kmplot_pfs(kmdf)}
              if(tolower(def_par) == "no" && tolower(plot_type) == "pfs"){graph <- kmplot_par_pfs(kmdf)}
              
              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            fit1 = surv_fit(f2, data = kmdf)
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv <- surv_pvalue(fit1)[[2]]
              if(pv <= p){
                png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              dev.off()
            }
          }
        } 
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }
  
  kmtable %>% filter(kmtable$`P-value` < p) %>% as.data.frame() -> sigKM
  list1 <- (list(kmtable, sigKM, kmdf))
  names(list1) <- c("kmtable", "sigKM", "transformeddf")
  sig_return <- readline(prompt = "Do you only want to save all the results, or just the significant results? (all/sig/none): ")
  if(tolower(sig_return) == "all" || tolower(sig_return) == "sig"){filen <- readline(prompt = "Please provide a name for your results: ")}
  
  if(tolower(sig_return) == "all"){
    res1 <- as.data.frame(list1[[1]])
    res2 <- as.data.frame(list1[[2]])
    name1 <- paste("./", filen, "_full_uniKM.xlsx", sep = "")
    name2 <- paste("./", filen, "_sig_uniKM.xlsx", sep = "")
    if(tolower(cutoff) == "optimal"){namec <- paste("./", filen, "_gene_cutoff.xlsx", sep = "")}

    write_xlsx(res1, name1)
    write_xlsx(res2, name2)
    if(tolower(cutoff) == "optimal"){write_xlsx(cuttable, namec)}
  }
  if(tolower(sig_return) == "sig"){
    res <- as.data.frame(list1[[2]])
    name1 <- paste("./", filen, "_sig_uniKM.xlsx", sep = "")
    if(tolower(cutoff) == "optimal"){namec <- paste("./", filen, "_gene_cutoff.xlsx", sep = "")}
    write_xlsx(res, name1)
    if(tolower(cutoff) == "optimal"){write_xlsx(cuttable, namec)}
  }
  return(list1)
  cat("All done :> \n")
}

multicox = function(kmdf){
  cat("Executing Multivariate Cox Analysis ... \n")
  message("For your query genes, do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
  input = readline(prompt = "Enter D or V (d/v): ")
  while(tolower(input) != "d" && tolower(input) != "v"){
    message("Answer not given, please try again :)")
    input = readline(prompt = "Enter D or V (d/v): ")
  }
  if(tolower(input) == "d"){
    message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
    genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
    genestostudy <- strsplit(genes, ", ")[[1]]
  }
  else if(tolower(input) == 'v'){
    message("Please provide a dataframe or a vector containing the gene names :)")
    gene_success <- FALSE
    while(!gene_success){
      genes <- readline(prompt = "Please enter the name of the dataframe or the vector of gene names: ")
      tryCatch({
        genestostudy <- get(genes)
        gene_success <- TRUE
      }, error = function(e){
        message("Error: Please provide a valid dataframe or a vector of gene names: ")
        gene_success <- FALSE
      })
    }
    
    cs_gts <- class(genestostudy)
    cs_gts <- cs_gts[1]
    if(cs_gts == "data.frame" || cs_gts == "data.table" || cs_gts == "tbl_df"){
      suppressWarnings({
        genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        while(is.na(genecolnum) == TRUE || genecolnum <= 0 || genecolnum > ncol(genestostudy)){
          message("Error detected, please try again :)")
          genecolnum <- as.numeric(readline(prompt = "Please provide the column number of the query genes: "))
        }
        colnames(genestostudy)[genecolnum] <- "geneid"
        genestostudy <- genestostudy$geneid
      })
    }
  }
  genestostudy <- gsub("-", "", genestostudy)
  kmdf %<>% as.data.frame()
  genestostudy <- genestostudy[genestostudy %in% colnames(kmdf)]
  
  type <- readline(prompt = "Are you looking into the overall survival or progression free survival? (os/pfs): ")
  cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
  while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal" && tolower(cutoff) != "continuous"){
    message("Answer not given, please try again :)")
    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
  }
  
  if(tolower(cutoff) == "median"){
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(genestostudy)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }
    
    for (i in unique(1:length(genestostudy))){
      z = genestostudy[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>% dplyr::select(all_of(l)) %>% unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      mdn = median(c)
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          c >= mdn,
          'High',
          'Low'
        )))} 
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    for (o in genestostudy){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
    }
    
    if(type == "os"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
    cox = coxph(f1, data = kmdf)
    zph = cox.zph(cox)
  }
  
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom"){
    if(tolower(cutoff) == "custom"){
      qval_low <- as.numeric(readline(prompt = "Please enter the low % cutoff between 0 and 1: "))
      qval_high <- as.numeric(readline(prompt = "Please enter the high % cutoff between 0 and 1: "))
    }
    
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(genestostudy)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }
    
    for (i in unique(1:length(genestostudy))){
      z = genestostudy[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>% dplyr::select(all_of(l)) %>% unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      if(tolower(cutoff) == "quartile"){quartile <- quantile(c, probs = c(0.25, 0.75))}
      if(tolower(cutoff) == "custom"){quartile <- quantile(c, probs = c(qval_low, qval_high))}
      
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          c >= quartile[[2]],
          'High',
          ifelse(
            c <= quartile[[1]], 
            "Low",
            "Medium"
          )
        )))
    }
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    for (o in genestostudy){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Medium')
    }
    
    if(type == "os"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
    cox = coxph(f1, data = kmdf)
    zph = cox.zph(cox)
  }
  
  if(tolower(cutoff) == "zscore"){
    zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(genestostudy)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }
    
    for (i in unique(1:length(genestostudy))){
      z = genestostudy[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>% dplyr::select(all_of(l)) %>% unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      zsc = scale(c)
      
      kmdf %<>% mutate("{name}" := as.factor(
        ifelse(
          zsc >= zcut,
          'High',
          ifelse(
            zsc <= -abs(zcut), 
            "Low",
            "Medium"
          )
        )))
    }
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    for (o in genestostudy){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Medium')
    }
    
    if(type == "os"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
    cox = coxph(f1, data = kmdf)
    zph = cox.zph(cox)
  }
  
  if(tolower(cutoff) == "optimal"){
    
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(genestostudy)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }
    
    if(type == "os"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    for (i in unique(1:length(genestostudy))){
      z = genestostudy[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>% dplyr::select(all_of(l)) %>% unlist() -> c
        fc <- as.formula(paste('survdf ~', paste(l)))
        tryCatch({
          mscutoff <- maxstat.test(fc, data = kmdf, smethod = "LogRank")
          o_cutoff <- mscutoff$estimate %>% as.double()
          
          kmdf %<>% mutate("{name}" := as.factor(
            ifelse(
              c >= o_cutoff,
              'High',
              'Low'
            )))
        }, error = function(e) {
          # Handle the error
          # Print an error message or take any other necessary actions
          print(paste("Error:", z, "contains non-zero cell counts within the specified proportion range, removed from the analysis"))
          # You can also choose to continue the loop without any further action
        })
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
    }
    
    kmdf <- kmdf %>% dplyr::select(where(~any(. != "High")))
    for (o in genestostudy){
      kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
    }
    
    f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
    cox = coxph(f1, data = kmdf)
    zph = cox.zph(cox)
  }
  
  if(tolower(cutoff) == "continuous"){
    if(type == "os"){survdf <- Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
    cox = coxph(f1, data = kmdf) 
  }
  
  unicolumn = c('Gene', 'Hazard Ratio', '95% CI', 'P-value', 'PH Check')
  coxtable = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })
  suppressWarnings({
    civalue <- as.numeric(readline(prompt = "Please enter your desired confidence interval between 0 and 1: "))
    while(is.na(civalue) == TRUE || civalue <= 0 || civalue >= 1){
      message("The CI value provided should be within 0 and 1, please try again :)")
      civalue <- as.numeric(readline(prompt = "Please enter your desired confidence interval between 0 and 1: "))
    }
  })
  
  ucoef = summary(cox)$coef[,2]
  upvalue = summary(cox)$coef[,5]
  ci = confint(cox, level = civalue)
  ci <- exp(ci)
  ci <- round(ci, digits = 2)
  uci <- paste(ci[,1], '~', ci[,2], sep = " ")
  phcheck <- zph$table[1:(length(rownames(zph$table)) - 1), 3]
  
  if(tolower(cutoff) == "median" || tolower(cutoff) == "optimal" || tolower(cutoff) == "continuous"){df = data.frame(genestostudy, ucoef, uci, upvalue, phcheck)}
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){
    dfname = rownames(summary(cox)$coef)
    dfname <- gsub("High", " High", dfname)
    dfname <- gsub("Low", " Low", dfname)
    df = data.frame(dfname, ucoef, uci, upvalue, phcheck)
  }
  coxtable <- rbind(coxtable, df)
  coxtable[length(rownames(coxtable)) + 1, 5] <- paste("Global = ", zph$table[length(rownames(zph$table)), 3])
  colnames(coxtable) = unicolumn
  rownames(coxtable) = NULL
  coxtable %>% filter(coxtable$`P-value` <= p) %>% as.data.frame() -> sigcox
  
  graph <- ggforest(cox, data = kmdf, refLabel = "Reference") %>% print()

  filen <- readline(prompt = "Please provide a name for your results: ")
  list1 <- (list(coxtable, sigcox))
  names(list1) <- c("full_cox", "sig_cox")
  res1 <- as.data.frame(list1[[1]])
  name1 <- paste("./", filen, "_full_multicox.xlsx", sep = "")
  
  write_xlsx(res1, name1)
  
  name3 <- paste("./", filen, "_forest_plot.png", sep = "")
  forestpar <- readline(prompt = "Do you want to change the size of the forest plot? (yes/no): ")
  if(tolower(forestpar) == "no" || tolower(forestpar) == "n"){ggsave(file = name3, width = 15, height = 15, dpi = 500)}
  if(tolower(forestpar) == "yes" || tolower(forestpar) == "y"){
    nw <- as.numeric(readline(prompt = "Please enter the value of the width: "))
    nh <- as.numeric(readline(prompt = "Please enter the value of the height: "))
    ggsave(file = name3, width = nw, height = nh, dpi = 500)
  }

  cat("Analysis Complete :> Hope you get something :> \n")
  return(cox)
}
