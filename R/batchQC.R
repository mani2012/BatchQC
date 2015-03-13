require(sva)
require(limma)
#' Checks for presence of batch effect and reports whether the batch needs to be adjusted
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @return pca Principal Components Analysis object of the data
#' @export
#' @examples
#' batchQC(data.matrix)
batchQC <- function(data.matrix, batch, mod=NULL)  {
  batchqc_heatmap(data.matrix, batch=batch, mod=mod)
  pca <- batchqc_pca(data.matrix, batch=batch, mod=mod)
  batchtest(pca, batch=batch, mod=mod)
  return(pca)
}
