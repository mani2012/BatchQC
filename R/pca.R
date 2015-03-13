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
    print ("Need to implement this part")
    # do something here
  }
  pca <- prcomp(t(data.matrix), retx=T, center=T, scale=T)
  pc <- pca$x 
  plot(pc[,1], pc[,2],col=cc)
  return(pca)
}
