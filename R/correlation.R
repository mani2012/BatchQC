#' Produce correlation heatmap plot
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @return Correlation heatmap plot
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, basemean=10000, ggstep=50, bbstep=2000, ccstep=800, 
#'     basedisp=100, bdispstep=-10, swvar=1000, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' pdata <- data.frame(batch, condition)
#' modmatrix = model.matrix(~as.factor(condition), data=pdata)
#' batchqc_correlation(data.matrix, batch, mod=modmatrix)
batchqc_correlation <- function(data.matrix, batch, mod = NULL) {
    lcpms <- log2CPM(data.matrix)
    lcounts <- lcpms$y
    
    cormat <- cor(lcounts)
    shinyInput <- getShinyInput()
    if (is.null(shinyInput)) {
        shinyInput <- list(data = data.matrix, batch = batch)
    }
    shinyInput <- c(shinyInput, list(cormat = cormat))
    setShinyInput(shinyInput)
    
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
    hcbHeatmap(cormat, ColSideColors = cc, RowSideColors = cc)
}

#' Produce Median Correlation plot
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @return Median Correlation plot
#' @import matrixStats
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, basemean=10000, ggstep=50, bbstep=2000, ccstep=800, 
#'     basedisp=100, bdispstep=-10, swvar=1000, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' pdata <- data.frame(batch, condition)
#' modmatrix = model.matrix(~as.factor(condition), data=pdata)
#' batchqc_corscatter(data.matrix, batch, mod=modmatrix)
batchqc_corscatter <- function(data.matrix, batch, mod = NULL) {
    
    lcpms <- log2CPM(data.matrix)
    lcounts <- lcpms$y
    
    cormat <- cor(lcounts)
    med_cor <- matrixStats::rowMedians(cormat)
    suggested_cutoff <- quantile(med_cor, p = 0.25) - 1.5 * IQR(med_cor)
    
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
    
    ylim <- c(pmin(min(med_cor), suggested_cutoff), max(med_cor))
    plot(med_cor, xaxt = "n", ylim = ylim, ylab = "median pairwise correlation",
        xlab = "", pch = 19, cex = 1.4, col = cc)
    axis(side = 1, at = seq(along = med_cor), labels = seq(along = med_cor))
    abline(v = seq(along = med_cor), lty = 3, lwd = 0.8)
    abline(h = suggested_cutoff, lty = 2)
} 
