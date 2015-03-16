#' Produce heatmap plots for the given data
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchqc_heatmap(data.matrix, batch)
batchqc_heatmap <- function(data.matrix, batch, mod=NULL)  {
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
  heatmap(data.matrix, ColSideColors=cc)  
}
