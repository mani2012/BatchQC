require("MCMCpack")

#' Generate simulated count data with batch effects for ngenes
#' 
#' @param ngenes Number of genes to simulate
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @param basemean Base mean
#' @param ggstep Gene to Gene step variation
#' @param bbstep Batch to Batch step variation
#' @param ccstep Condition to Condition step variation
#' @param basedisp Base Dispersion
#' @param bdispstep Batch to Batch Dispersion step variation
#' @param swvar Sample-wise extra variation
#' @param seed Random seed for reproducibility
#' @return RNA Seq count data matrix
#' @import MCMCpack
#' @export
#' @examples
#' rnaseq_sim()
#' rnaseq_sim(ngenes=100, nbatch=5, seed=1234)
#' rnaseq_sim(ngenes=100, nbatch=3, ncond=2, npercond=10, basemean=10000,
#'     ggstep=50, bbstep=20000, ccstep=8000, basedisp=100, bdispstep=10, 
#'     swvar=1000, seed=1234)
rnaseq_sim <- function(ngenes = 50, nbatch = 3, ncond = 2, npercond = 10, 
    basemean = 10000, ggstep = 50, bbstep = 20000, ccstep = 8000, 
    basedisp = 100, bdispstep = 10, swvar = 1000, seed = 1000) {
    set.seed(seed, kind = NULL, normal.kind = NULL)
    mu <- seq(0, length.out = ngenes, by = ggstep)
    bmu <- seq(0, length.out = nbatch, by = bbstep)
    cmu <- seq(0, length.out = ncond, by = ccstep)
    bsize <- seq(basedisp, length.out = nbatch, by = bdispstep)
    ncol <- nbatch * ncond * npercond
    A.matrix <- matrix(0, nrow = ngenes, ncol = ncol)
    samplewisevar <- swvar*rbeta(ncol,2,2)
    for (i in 1:ngenes) {
        genei <- c()
        for (j in 1:nbatch) {
            for (k in 1:ncond) {
                for (l in 1:npercond) {
                    genei <- c(genei, rnbinom(1, size = bsize[j], 
                        mu = basemean+mu[i]+bmu[j]+cmu[k]+samplewisevar[j*k*l]))
                }
            }
        }
        A.matrix[i, ] <- genei
    }
    return(A.matrix)
} 
