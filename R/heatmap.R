#' Produce heatmap plots for the given data
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @export
#' @examples
#' batchqc_heatmap(data.matrix)
#' batchqc_heatmap(data.matrix, nbatch=5)
#' batchqc_heatmap(data.matrix, nbatch=3, ncond=2, npercond=10)
batchqc_heatmap <- function(data.matrix, nbatch=3, ncond=2, npercond=10)  {
  bc <- rainbow(nbatch)
  cc <- rep(bc, each=ncond*npercond)
  heatmap(data.matrix, ColSideColors=cc)  
}
