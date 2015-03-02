#' Performs test to check whether batch needs to be adjusted 
#' 
#' @param pca PCA object from principal component analysis
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @export
#' @examples
#' batchtest(pca)
#' batchtest(pca, nbatch=5)
#' batchtest(pca, nbatch=3, ncond=2, npercond=10)
batchtest <- function(pca, nbatch=3, ncond=2, npercond=10)  {
  pc <- pca$x 
  batch <- rep(1:nbatch, each=ncond*npercond)
  summary(lm(pc[,1]~batch))
  summary(lm(pc[,2]~batch))
  summary(lm(pc[,3]~batch))
  summary(lm(pc[,4]~batch))
  summary(lm(pc[,5]~batch))
}
