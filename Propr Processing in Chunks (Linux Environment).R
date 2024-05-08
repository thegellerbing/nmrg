propr.chunk <- function(counts, metric = c("rho", "phi", "phs", "cor", "vlr"),
                        ivar = NA, symmetrize = FALSE, alpha=NA, p=100, fdr=0.05,
                        n=100, ncores = 1, interval=seq(0 , 1, 0.05),
                        dir = NA){
  
  # divide data into chunks
  # I recommend using large chunk size n,
  # otherwise dividing the data into too many chunks will slow down the computation
  l <- blocks2combs(counts, n)
  combs <- l[[1]]
  split <- l[[2]]
  
  # parallelize propr computation for i chunks
  doMC::registerDoMC(cores = ncores)
  `%dopar%` <- foreach::`%dopar%`
  
  RES <- foreach::foreach(i = 1:nrow(combs)) %dopar% {
    
    # get chunk
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    chunk = subset(counts, select = c(batch1, batch2))
    
    # compute propr 
    l_propr <- chunk2propr(i, chunk, metric=metric, ivar=ivar, symmetrize=symmetrize, alpha=alpha, p=p, interval=interval, fdr=fdr)
    
    if(is.na(dir)){
      # return propr matrix, fdr
      l_propr
    }else{
      # save data if required
      file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
      print(paste("----saving tmp file[", i, "] to ", file2, sep=""))
      write.csv(l_propr[[1]], file=file2)
      
      # return cutoff and FDR
      list(NULL, l_propr[[2]])
    }
  }
  
  # collect files
  if(!is.na(dir)){
    RES <- file2res(RES, combs, dir)
  }
  
  # merge chunks
  l_propr <- chunk2full(counts, RES, split, combs, fdr)
  
  return(l_propr)
}

blocks2combs <- function(counts, n){
  
  # define blocks
  # the resulting chunks will have size of at most [n, n]
  nblocks = ncol(counts) %/% n
  ngroup <- rep(1:nblocks, each = n)
  leftover <- ncol(counts) - length(ngroup)
  if(leftover > 0) ngroup <- c(ngroup, rep(nblocks + 1, leftover))
  
  # check size
  # here I decided to stop computation if no more than 2 groups are generated (so only 1 chunk)
  # so if this happens you better change n
  if (length(unique(ngroup)) <= 2){
    stop(paste("ERROR: chunk size ", n, " is too big for data frame of size [", nrow(counts), "][", ncol(counts), "]", sep=""))
  }
  
  # split groups
  # each row in combs define a chunk
  split <- split(1:ncol(counts), ngroup)
  combs <- expand.grid(1:length(split), 1:length(split))
  combs <- t(apply(combs, 1, sort))
  combs <- unique(combs) 
  combs <- combs[combs[,1] != combs[,2],]
  
  print(paste("----runing propr for ", nrow(combs), " chunks of size ", n))
  
  return(list(combs, split))
}

chunk2propr <- function(i, chunk, metric="rho", ivar=NA, symmetrize=FALSE,
                        alpha, p=100, interval=seq(0 , 1, 0.05), fdr=0.05 ){
  
  print(i)
  
  # compute propr for chunk
  rho.i <- propr(chunk, metric = metric, ivar = ivar, alpha = alpha, p=p)
  rho.i <- updateCutoffs(rho.i,interval)
  
  # get cutoff | fdr
  df <- data.frame('cutoff'=rho.i@fdr[,'cutoff'], 'FDR'=rho.i@fdr[,'FDR'])
  
  return(list(rho.i@matrix, df))
}

chunk2full <- function(counts, RES, split, combs, fdr){
  
  # define variables
  QUILT <- matrix(0, ncol(counts), ncol(counts))
  d <- data.frame('cutoff'=RES[[1]][[2]][,'cutoff'])
  
  # collect chunks
  for(i in 1:nrow(combs)){
    
    # add fdr.i
    d[paste('FDR', i, sep="")] <- RES[[i]][[2]][,'FDR']
    
    # Fill final matrix with each chunk
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    patch.i <- c(batch1, batch2)
    QUILT[patch.i, patch.i] <- RES[[i]][[1]]
  }
  
  # average fdr
  df <- data.frame('cutoff'=RES[[1]][[2]][,'cutoff'], 'FDR'=round(rowSums(d[,2:ncol(d)])/i,4))
  # cutoff
  cutoff <- min(df[df[,"FDR"]<fdr,"cutoff"])
  
  # rename columns & rows
  matrix <- QUILT
  rownames(matrix) <- colnames(counts)
  colnames(matrix) <- colnames(counts)
  
  return(list(matrix, cutoff, df))
}

file2res <- function(RES, combs, dir){
  
  print(paste("----collecting ", nrow(combs), " chunk files previously saved in ", dir, sep=""))
  
  for(i in 1:nrow(combs)){
    file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
    csv <- read.csv(file2, row.names = 1, header= TRUE)
    matrix <- data.matrix(csv)
    RES[[i]][[1]] <- matrix
  }
  
  return(RES)
}
