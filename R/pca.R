#' Performs principal component analysis and 
#' produces plot of the first two principal components
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @return PCA object from principal component analysis
#' @export
#' @examples
#' batchqc_pca(data.matrix)
#' batchqc_pca(data.matrix, nbatch=5)
#' batchqc_pca(data.matrix, nbatch=3, ncond=2, npercond=10)
batchqc_pca <- function(data.matrix, nbatch=3, ncond=2, npercond=10)  {
  bc <- rainbow(nbatch)
  cc <- rep(bc, each=ncond*npercond)
  pca <- prcomp(t(data.matrix), retx=T, center=T, scale=T)
  pc <- pca$x 
  plot(pc[,1], pc[,2],col=cc)
  return(pca)
}
