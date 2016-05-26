#' Performs test to check whether batch needs to be adjusted 
#' 
#' @param pca PCA object from principal component analysis
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @return Summary of linear regression of first five principal components 
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
#' pca <- batchqc_pca(data.matrix, batch, mod=modmatrix)
#' batchtest(pca, batch, mod=modmatrix)
batchtest <- function(pca, batch, mod = NULL) {
    pc <- pca$x
    fbatch <- as.factor(batch)
    condition <- NULL
    if (!is.null(mod)) {
        condition = mod[, 2]
    }
    if (is.null(condition)) {
        s1 <- summary(lm(pc[, 1] ~ fbatch))
        s2 <- summary(lm(pc[, 2] ~ fbatch))
        s3 <- summary(lm(pc[, 3] ~ fbatch))
        s4 <- summary(lm(pc[, 4] ~ fbatch))
        s5 <- summary(lm(pc[, 5] ~ fbatch))
    } else {
        fcond <- as.factor(condition)
        s1 <- summary(lm(pc[, 1] ~ fbatch + fcond))
        s2 <- summary(lm(pc[, 2] ~ fbatch + fcond))
        s3 <- summary(lm(pc[, 3] ~ fbatch + fcond))
        s4 <- summary(lm(pc[, 4] ~ fbatch + fcond))
        s5 <- summary(lm(pc[, 5] ~ fbatch + fcond))
    }
    print(s1)
    print(s2)
    print(s3)
    print(s4)
    print(s5)
    retval <- c(s1, s2, s3, s4, s5)
    return(retval)
} 
