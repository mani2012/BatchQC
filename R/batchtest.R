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
  s1 <- summary(lm(pc[,1]~batch))
  s2 <- summary(lm(pc[,2]~batch))
  s3 <- summary(lm(pc[,3]~batch))
  s4 <- summary(lm(pc[,4]~batch))
  s5 <- summary(lm(pc[,5]~batch))
  print(s1)
  print(s2)
  print(s3)
  print(s4)
  print(s5)
}
