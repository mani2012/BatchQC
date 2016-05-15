#' Produce heatmap plots for the given data
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @param max_display Maximum number of rows to display in heat map 
#' @return Heatmap plots for the given data 
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' pdata <- data.frame(batch, condition)
#' modmatrix = model.matrix(~as.factor(condition), data=pdata)
#' batchqc_heatmap(data.matrix, batch, mod=modmatrix)
batchqc_heatmap <- function(data.matrix, batch, mod = NULL, max_display = 50) {
    size <- dim(data.matrix)[1]
    reduced.data.matrix <- data.matrix
    if (size > max_display) {
        reduced.data.matrix <- data.matrix[1:max_display, ]
    }
    lcpm <- log2CPM(reduced.data.matrix)
    lcounts <- lcpm$y
    
    fbatch <- as.factor(batch)
    nbatch <- nlevels(fbatch)
    bc <- rainbow(nbatch)
    intbatch <- as.integer(fbatch)
    colorfun <- function(i) {
        return(bc[i])
    }
    cc <- sapply(intbatch, colorfun, simplify = TRUE)
    if (!is.null(mod)) {
        # print ('Need to implement this part') do something here
    }
    hcbHeatmap(lcounts, ColSideColors = cc)
} 
