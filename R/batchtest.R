#' Performs test to check whether batch needs to be adjusted 
#' 
#' @param pca PCA object from principal component analysis
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchtest(pca, batch)
batchtest <- function(pca, batch, mod=NULL)  {
  pc <- pca$x 
  summary(lm(pc[,1]~batch))
  summary(lm(pc[,2]~batch))
  summary(lm(pc[,3]~batch))
  summary(lm(pc[,4]~batch))
  summary(lm(pc[,5]~batch))
}
