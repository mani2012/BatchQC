require("MCMCpack")

#' Generate simulated count data with batch effects for ngenes
#' 
#' @param ngenes Number of genes to simulate
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @param ggstep Gene to Gene step variation
#' @param bbstep Batch to Batch step variation
#' @param ccstep Condition to Condition step variation
#' @param bvarstep Batch to Batch variance step variation
#' @param seed Random seed for reproducibility
#' @return RNA Seq count data matrix
#' @export
#' @examples
#' rnaseq_sim()
#' rnaseq_sim(ngenes=100, nbatch=5, seed=1234)
#' rnaseq_sim(ngenes=100, nbatch=3, ncond=2, npercond=10, ggstep=5, 
#'            bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
rnaseq_sim <- function(ngenes=50, nbatch=3, ncond=2, npercond=10, ggstep=5, 
                       bbstep=15000, ccstep=10000, bvarstep=2, seed=1000)  {
  set.seed(seed, kind = NULL, normal.kind = NULL)
  mu <- seq(ggstep, length.out=ngenes, by=ggstep)
  bmu <- seq(bbstep, length.out=nbatch, by=bbstep)
  cmu <- seq(ccstep, length.out=ncond, by=ccstep)
  bsize <- seq(bvarstep, length.out=nbatch, by=bvarstep)
  A.matrix <- matrix(0, nrow=ngenes, ncol=nbatch*ncond*npercond)
  for (i in 1:ngenes)  {
    genei <- c()
    for (j in 1:nbatch)  {
      batchmu <- rnorm(1, mean=bmu[j], sd=1)
      size <- rinvgamma(1, shape=bsize[j], scale=1)
      for (k in 1:ncond)  {
        condmu <- rnorm(1, mean=cmu[k], sd=1)
        genei <- c(genei, rnbinom(npercond,size=size,mu=mu[i]+batchmu+condmu))
      }
    }
    A.matrix[i,] <- genei
  }
  return(A.matrix)
}
