require(sva)
require(limma)
#' Checks for presence of batch effect and adjusts for it
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @return clean.data Cleaned data
#' @export
#' @examples
#' batchqc_pipeline(data.matrix)
#' batchqc_pipeline(data.matrix, nbatch=3, ncond=2, npercond=10)
batchqc_pipeline <- function(data.matrix, nbatch=3, ncond=2, npercond=10)  {
  batchqc_heatmap(data.matrix, nbatch=3, ncond=2, npercond=10)
  pca <- batchqc_pca(data.matrix, nbatch=3, ncond=2, npercond=10)
  batchtest(pca, nbatch=3, ncond=2, npercond=10)

  ### apply combat
  batch <- rep(1:nbatch, each=ncond*npercond)
  sample <- 1:nbatch*ncond*npercond
  pdata <- data.frame(sample, batch)
  modcombat = model.matrix(~1, data=pdata)  
  combat_data.matrix = ComBat(dat=data.matrix, batch=batch, mod=modcombat)

  batchqc_heatmap(combat_data.matrix, nbatch=3, ncond=2, npercond=10)
  combat_pca <- batchqc_pca(combat_data.matrix, nbatch=3, ncond=2, npercond=10)
  batchtest(combat_pca, nbatch=3, ncond=2, npercond=10)
  
  return(combat_data.matrix)
}
data.matrix <- rnaseq_sim(ngenes=50, nbatch=3, ncond=2, npercond=10, ggstep=5,
                          bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
clean.data <- batchqc_pipeline(data.matrix, nbatch=3, ncond=2, npercond=10)

