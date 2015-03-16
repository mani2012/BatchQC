#' Performs principal component analysis and 
#' produces plot of the first two principal components
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @return PCA object from principal component analysis
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchqc_pca(data.matrix, batch)
batchqc_pca <- function(data.matrix, batch, mod=NULL)  {
  fbatch <- as.factor(batch)
  nbatch <- nlevels(fbatch)
  bc <- rainbow(nbatch)
  intbatch <- as.integer(fbatch)
  colorfun <- function(i) { return(bc[i]) }
  cc <- sapply(intbatch, colorfun, simplify=TRUE)
  if (!is.null(mod)) {
    #print ("Need to implement this part")
    # do something here
  }
  pca <- prcomp(t(data.matrix), retx=T, center=T, scale=T)
  pc <- pca$x 
  xlab <- "Principal Component PC1"
  ylab <- "Principal Component PC2"
  plot(pc[,1], pc[,2],col=cc, xlab=xlab, ylab=ylab)
  return(pca)
}

#' Performs PCA svd variance decomposition and 
#' produces plot of the first two principal components
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @return res PCA list with two components v and d.
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchqc_pca_svd(data.matrix, batch)
batchqc_pca_svd <- function(data.matrix, batch, mod=NULL)  {
  res <- makeSVD(data.matrix)
  pcRes(res$v,res$d, modmatrix[,2], batch)
  fcond <- as.factor(modmatrix[,2])
  ncond <- nlevels(fcond)
  bc <- rainbow(ncond)
  intcond <- as.integer(fcond)
  colorfun <- function(i) { return(bc[i]) }
  cc <- sapply(intcond, colorfun, simplify=TRUE)
  #fbatch <- as.factor(batch)
  #nbatch <- nlevels(fbatch)
  #bc <- rainbow(nbatch)
  #intbatch <- as.integer(fbatch)
  #colorfun <- function(i) { return(bc[i]) }
  #cc <- sapply(intbatch, colorfun, simplify=TRUE)
  plotPC(res$v,res$d, 
         col=cc, # color by condition
         pch=19, main="PCA plot",
         xlim=c(min(res$v[,1])-.08,max(res$v[,1])+.08),
         ylim=c(min(res$v[,2])-.08,max(res$v[,2])+.08))
  text(res$v[,1], res$v[,2], batch, pos=1, cex=0.6)
  return(res)
}

#' Compute singular value decomposition
#'
#' @param x matrix of genes by sample (ie. the usual data matrix)
#' @return returns a list of svd components v and d
#' @export
makeSVD <- function(x){
  x <- as.matrix(x)
  s <- fast.svd(x-rowMeans(x))
  
  v <- s$v
  rownames(v) <- colnames(x)
  
  s <- list(v=v, d=s$d)
  return(s)
}


#' Compute variance of each principal component and how they correlate with batch and cond  
#'
#' @param v from makeSVD
#' @param d from makeSVD
#' @param condition factor describing experiment
#' @param batch factor describing batch
#' @return A dataframe containig variance, cum. variance, cond.R-sqrd, batch.R-sqrd
#' @export
pcRes <- function(v, d, condition=NULL, batch=NULL){
  pcVar <- round((d^2)/sum(d^2)*100,2)
  cumPcVar <- cumsum(pcVar) 
  
  if(!is.null(condition)){
    cond.R2 <- function(y) round(summary(lm(y~condition))$r.squared*100,2)
    cond.R2 <- apply(v, 2, cond.R2) 
  }
  
  if(!is.null(batch)){
    batch.R2 <- function(y) round(summary(lm(y~batch))$r.squared*100,2)
    batch.R2 <- apply(v, 2, batch.R2)
  }
  
  if(is.null(condition) & is.null(batch)){
    res <- data.frame(propVar=pcVar, 
                      cumPropVar=cumPcVar)
  }
  
  if(!is.null(batch) & is.null(condition)){
    res <- data.frame(propVar=pcVar, 
                      cumPropVar=cumPcVar, 
                      batch.R2=batch.R2)
  }
  
  if(!is.null(condition) & is.null(batch)){
    res <- data.frame(propVar=pcVar, 
                      cumPropVar=cumPcVar, 
                      cond.R2=cond.R2)
  }
  
  if(!is.null(condition) & !is.null(batch)){
    res <- data.frame(propVar=pcVar, 
                      cumPropVar=cumPcVar, 
                      cond.R2=cond.R2,
                      batch.R2=batch.R2)
  }
  print(res)
  return(res)
}

#' Plot first 2 principal components
#'
#' @param v from makeSVD
#' @param d from makeSVD
#' @param ... pass options to internal plot fct.
#' @return a plot 
#' @export
plotPC <- function(v, d, ...){
  pcVar <- round((d^2)/sum(d^2)*100,2)
  
  xl <- sprintf("pc 1: %.2f%% variance", pcVar[1])  
  yl <- sprintf("pc 2: %.2f%% variance", pcVar[2]) 
  
  plot(v[,1], v[,2], xlab=xl, ylab=yl, ...)
}
