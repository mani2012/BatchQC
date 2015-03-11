require(sva)
require(limma)
#' Checks for presence of batch effect and reports whether the batch needs to be adjusted
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param nbatch Number of batches in the data
#' @param ncond Number of conditions in the data
#' @param npercond Number of samples per condition per batch in the data
#' @return pca Principal Components Analysis object of the data
#' @export
#' @examples
#' batchQC(data.matrix)
#' batchQC(data.matrix, nbatch=3, ncond=2, npercond=10)
batchQC <- function(data.matrix, nbatch=3, ncond=2, npercond=10)  {
  batchqc_heatmap(data.matrix, nbatch=3, ncond=2, npercond=10)
  pca <- batchqc_pca(data.matrix, nbatch=3, ncond=2, npercond=10)
  batchtest(pca, nbatch=3, ncond=2, npercond=10)
  return(pca)
}
