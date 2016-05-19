library(BatchQC)
context("BatchQC functionality")

test_that("batchQC_analyze", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    pca <- batchQC_analyze(data.matrix, batch)
    expect_equal(length(pca$x), 3000)
})

test_that("batchQC", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond,
        npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2,
        seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    outputfile <- batchQC(data.matrix, batch=batch, condition=condition, 
        view_report=FALSE, interactive=FALSE)
    expect_equal(basename(outputfile), "batchqc_report.html")
})

test_that("combatPlot", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond,
        npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2,
        seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    mod = model.matrix(~as.factor(condition), data = pdata)
    kstest <- combatPlot(data.matrix, batch, mod=mod)
    expect_equal(signif(kstest$p.value, 4), 0.3187)
})

test_that("batchQC_condition_adjusted", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    condition_adjusted <- batchQC_condition_adjusted(data.matrix, batch, 
        condition)
    expect_equal(dim(data.matrix), dim(condition_adjusted))
})

test_that("batchqc_explained_variation", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    batchqc_ev <- batchqc_explained_variation(data.matrix, condition, batch)
    expect_equal(dim(batchqc_ev$explained_variation), c(50,3))
})

test_that("batchqc_pc_explained_variation", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    pca <- batchqc_pca(data.matrix, batch, mod=modmatrix)
    pcs <- t(data.frame(pca$x))
    explained_variation <- batchqc_pc_explained_variation(pcs, pca$sdev^2, 
        condition, batch)
    expect_equal(dim(explained_variation), c(50,7))
})

test_that("log2CPM", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond,
        npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2,
        seed=1234)
    data.matrix <- as.matrix(data.matrix)
    ld <- log2CPM(data.matrix)
    expect_equal(dim(ld$y), dim(data.matrix))
    expect_equal(length(ld$lib.size), dim(data.matrix)[2])
})

test_that("batchQC_num.sv", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    num.sv <- batchQC_num.sv(data.matrix, modmatrix)
    expect_equal(num.sv, 2)
})

test_that("batchQC_fsva_adjusted", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    sva.object <- batchQC_sva(data.matrix, mod=modmatrix)
    svaf_data <- batchQC_fsva_adjusted(data.matrix, modmatrix, sva.object)
    expect_equal(dim(svaf_data), dim(data.matrix))
})

test_that("batchQC_svregress_adjusted", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    sva.object <- batchQC_sva(data.matrix, mod=modmatrix)
    svar_data <- batchQC_svregress_adjusted(data.matrix, modmatrix, sva.object)
    expect_equal(dim(svar_data), dim(data.matrix))
})

test_that("batchQC_shapeVariation", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond,
        npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2,
        seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    sl <- batchQC_shapeVariation(data.matrix, groups=batch)
    expect_equal(length(sl), 8)
})

test_that("batchqc_pca", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    pca <- batchqc_pca(data.matrix, batch, mod=modmatrix)
    expect_equal(length(pca$x), 3000)
})

test_that("batchqc_pca_svd", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    res1 <- batchqc_pca_svd(data.matrix, batch, mod=modmatrix)
    expect_equal(dim(res1), c(50,4))
})

test_that("batchtest", {
    nbatch <- 3
    ncond <- 2
    npercond <- 10
    data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
        npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    batch <- rep(1:nbatch, each=ncond*npercond)
    condition <- rep(rep(1:ncond, each=npercond), nbatch)
    pdata <- data.frame(batch, condition)
    modmatrix = model.matrix(~as.factor(condition), data=pdata)
    pca <- batchqc_pca(data.matrix, batch, mod=modmatrix)
    retval <- batchtest(pca, batch, mod=modmatrix)
    expect_equal(length(retval), 55)
})

test_that("rnaseq_sim", {
    data.matrix <- rnaseq_sim(ngenes=100, nbatch=3, ncond=2, npercond=10, 
        ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
    expect_equal(dim(data.matrix), c(100,60))
})

