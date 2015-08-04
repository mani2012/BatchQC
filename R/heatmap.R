#' Produce heatmap plots for the given data
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param max_display Maximum number of rows to display in heat map 
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchqc_heatmap(data.matrix, batch)
batchqc_heatmap <- function(data.matrix, batch, mod=NULL, max_display=50)  {
  size <- dim(data.matrix)[1]
  reduced.data.matrix <- data.matrix
  if (size > max_display)  {
    reduced.data.matrix <- data.matrix[1:max_display,]
  }
  lcpm <- log2CPM(reduced.data.matrix)
  lcounts <- lcpm$y
  
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
  hcbHeatmap(lcounts, ColSideColors=cc)  
}
