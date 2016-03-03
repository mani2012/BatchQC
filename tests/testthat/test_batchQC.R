library(BatchQC)
context("BatchQC functionality")

test_that("BatchQC functionality", {
  nbatch <- 3
  ncond <- 2
  npercond <- 10
  data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=npercond, 
                            ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
  batch <- rep(1:nbatch, each=ncond*npercond)
  pca <- batchQC_analyze(data.matrix, batch)
  expect_equal(length(pca$x), 3000)
})

