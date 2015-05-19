#' Produce correlation heatmap plot
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchqc_correlation(data.matrix, batch)
batchqc_correlation <- function(data.matrix, batch, mod=NULL) {
  lcpms <- log2CPM(data.matrix)
  lcounts <- lcpms$y
  
  cormat <- cor(lcounts)
  
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
  hcbHeatmap(cormat, ColSideColors=cc, RowSideColors=cc)
}

batchqc_corscatter <- function(data.matrix, batch, mod=NULL) {
  
  lcpms <- log2CPM(data.matrix)
  lcounts <- lcpms$y
  
  cormat <- cor(lcounts)
  med_correlation <- matrixStats::rowMedians(cormat)
  
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
  
}
