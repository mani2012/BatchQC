# Functions to perform surrogate variable analysis within the
# BatchQC Package

require(sva)


###### Estimate Surrogate Variables ######

#' Returns the number of surrogate variables to estimate in the model using 
#' a permutation based procedure
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param modmatrix Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @return Number of Surrogate variables found
#' @import sva
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
#' batchQC_num.sv(data.matrix, modmatrix)
batchQC_num.sv = function(data.matrix, modmatrix) {
    data.matrix = as.matrix(data.matrix)
    modmatrix = as.matrix(modmatrix)
    num.sv.be = sva::num.sv(dat = data.matrix, mod = modmatrix, 
        method = "be", vfilter = NULL, B = 20, seed = 47)
    return(num.sv.be)
}

## Output the data with (known) batch effect removed but
## biological heterogeneity preserved
batchQC_psva = function(data.matrix, batch) {
    data.matrix = as.matrix(data.matrix)
    batch = as.factor(batch)
    psva.output = psva(dat = data.matrix, batch = batch)
    return(psva.output)
}

#' Estimate the surrogate variables using the 2 step approach proposed by 
#' Leek and Storey 2007 
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param modmatrix Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @return Surrogate variables analysis object
#' @import sva
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
#' batchQC_sva(data.matrix, modmatrix)
batchQC_sva = function(data.matrix, modmatrix) {
    n.sv = batchQC_num.sv(data.matrix, modmatrix)
    modmatrix0 = model.matrix(~1, data = data.frame(t(data.matrix)))
    sva.object = sva::sva(dat = data.matrix, mod = modmatrix, 
        mod0 = modmatrix0, n.sv = n.sv)
    return(sva.object)
}


###### Clean Data by Removing Surrogate Variables ######

#' Use frozen surrogate variable analysis to remove the surrogate variables 
#' inferred from sva 
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param modmatrix Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @param sva.object SVA object
#' @return Frozen Surrogate variables adjusted data
#' @import sva
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
#' sva.object <- batchQC_sva(data.matrix, mod=modmatrix)
#' batchQC_fsva_adjusted(data.matrix, modmatrix, sva.object)
batchQC_fsva_adjusted = function(data.matrix, modmatrix, sva.object) {
    fsva.object = sva::fsva(data.matrix, modmatrix, sva.object, 
        data.matrix, method = "fast")
    fsva.adjusted.db = fsva.object$db
    return(fsva.adjusted.db)
}

#' Regress the surrogate variables out of the expression data 
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param modmatrix Model matrix for outcome of interest and other covariates 
#'     besides batch
#' @param sva.object SVA object
#' @return Surrogate variables regress adjusted data
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
#' sva.object <- batchQC_sva(data.matrix, mod=modmatrix)
#' batchQC_svregress_adjusted(data.matrix, modmatrix, sva.object)
batchQC_svregress_adjusted = function(data.matrix, modmatrix, 
    sva.object) {
    y = data.matrix
    X = cbind(modmatrix, sva.object$sv)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    P = ncol(modmatrix)
    svregress.adjusted.db = y - t(as.matrix(X[, -c(1:P)]) %*% 
        beta[-c(1:P), ])
    return(svregress.adjusted.db)
} 
