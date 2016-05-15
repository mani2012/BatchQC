#' Perform Mean and Variance batch variation analysis
#' 
#' @param data Given data
#' @param groups a character vector indicating sample group membership
#' @param plot Indicate whether to generate plot
#' @param groupCol group color
#' @return Mean and Variance batch variation Overall and Pairwise p-values
#' @import gplots moments
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, 
#'     npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, 
#'     seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' batchQC_shapeVariation(data.matrix, groups=batch)
batchQC_shapeVariation = function(data, groups, plot = FALSE, 
    groupCol = NULL) {
    
    groupsorder <- order(groups)
    sortdata <- data[,groupsorder]
    sortgroups <- groups[groupsorder]
    groupCol <- groupCol[groupsorder]
    gnormdata <- gnormalize(sortdata)
    #Y = fitMoments(sortdata)
    Y = fitMoments(gnormdata)
    
    ps <- overallPvalue(Y, sortgroups)
    #batch_ps <- batchEffectPvalue(sortdata, sortgroups)
    batch_ps <- batchEffectPvalue(gnormdata, sortgroups)
    
    mpval = ps[1]
    vpval = ps[2]
    spval = ps[3]
    kpval = ps[4]
    mpvaltext <- round(mpval, 4)
    if (mpvaltext == 0) {
        mpvaltext <- "< 0.0001"
    }
    vpvaltext <- round(vpval, 4)
    if (vpvaltext == 0) {
        vpvaltext <- "< 0.0001"
    }
    spvaltext <- round(spval, 4)
    if (spvaltext == 0) {
        spvaltext <- "< 0.0001"
    }
    kpvaltext <- round(kpval, 4)
    if (kpvaltext == 0) {
        kpvaltext <- "< 0.0001"
    }
    
    mpval2 = batch_ps[1]
    vpval2 = batch_ps[2]
    spval2 = batch_ps[3]
    kpval2 = batch_ps[4]
    mpvaltext2 <- round(mpval2, 4)
    if (mpvaltext2 == 0) {
        mpvaltext2 <- "< 0.0001"
    }
    vpvaltext2 <- round(vpval2, 4)
    if (vpvaltext2 == 0) {
        vpvaltext2 <- "< 0.0001"
    }
    spvaltext2 <- round(spval2, 4)
    if (spvaltext2 == 0) {
        spvaltext2 <- "< 0.0001"
    }
    kpvaltext2 <- round(kpval2, 4)
    if (kpvaltext2 == 0) {
        kpvaltext2 <- "< 0.0001"
    }
    
    if (plot) {
        
        pal = colorRampPalette(c("red", "orange", "white", "steelblue3", 
            "navy"))(n = 99)
        
        main = paste0(
            "\n\nBatch Effect Variation Analysis \n Mean p-value: Overall = ", 
            mpvaltext, ", Pairwise = ", mpvaltext2, 
            "\n Variance p-value: Overall = ", vpvaltext, ", Pairwise = ", 
            vpvaltext2, "\n Skewness p-value: Overall = ", spvaltext, 
            ", Pairwise = ", spvaltext2, "\n Kurtosis p-value: Overall = ", 
            kpvaltext, ", Pairwise = ", kpvaltext2)
        
        gplots::heatmap.2(t(Y), trace = "none", Rowv = FALSE, 
            Colv = FALSE, dendrogram = "none", col = pal, ColSideColors = 
            groupCol, density.info = "none", scale = "row", cexRow = 1.15, 
            colsep = cumsum(table(sortgroups)), main = main)
        
        legend("bottomleft", legend = unique(sortgroups), pch = 19, 
            col = unique(groupCol), title = "Batch")
    }
    c(ps, batch_ps)
}

fitMoments <- function(x) {
    # Compute Mean and Variance
    smean <- apply(x, 2, mean)
    svariance <- apply(x, 2, var)
    sskew <- apply(x, 2, skewness)
    skurt <- apply(x, 2, kurtosis)
    sfit <- cbind(smean, svariance, sskew, skurt)
    colnames(sfit) = c("Mean", "Variance", "Skew", "Kurtosis")
    rownames(sfit) = colnames(x)
    return(sfit)
}

#' Perform Genewise Normalization of the given data matrix
#' 
#' @param dat Given data matrix
#' @return gnormdata Genewise Normalized data matrix
#' @export
#' @examples
#' dat <- matrix(1:10, 2)
#' gnormdata <- gnormalize(dat)
gnormalize <- function(dat) {
    # Compute Genewise Mean and Standard deviation
    smean <- apply(dat, 1, mean)
    ssd <- apply(dat, 1, sd)
    smeanmat <- apply(row(dat), c(1,2), function(x, smean) smean[x], smean)
    ssdmat <- apply(row(dat), c(1,2), function(x, ssd) ssd[x], ssd)
    # Genewise normalize the data
    gnormdata <- matrix(mapply(
        function(x, smean, ssd) ((x-smean)/ssd) + smean, 
        dat, smeanmat, ssdmat), nrow = nrow(dat))
    return(gnormdata)
}

overallPvalue <- function(Y, groups) {
    mod1 = model.matrix(~as.factor(groups), data = data.frame(Y))
    # model with batch
    mod0 = model.matrix(~1, data = data.frame(Y))  # reduced model 
    
    ### Standard F test ###
    n <- dim(Y)[1]
    Id <- diag(n)
    df1 <- dim(mod1)[2]
    df0 <- dim(mod0)[2]
    
    resid <- t(Y) %*% (Id - mod1 %*% solve(t(mod1) %*% mod1) %*% 
        t(mod1))
    # residuals full model
    rss1 <- rowSums(resid * resid)  ## SSE full model
    
    resid0 <- t(Y) %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% 
        t(mod0))  # residuals reduced model
    rss0 <- rowSums(resid0 * resid0)  ## SSE reduced model
    
    #delta <- (apply(t(Y), 1, mean) * 0.01)^2
    delta <- 0
    Fstat = ((rss0 - rss1)/(df1 - df0))/(delta + rss1/(n - df1))
    # Compute p-value
    ps = 1 - pf(Fstat, df1 - df0, n - df1)
    return(ps)
}

delta_f.pvalue <- function(dat, mod, mod0) {
    ## F-test (full/reduced model) and returns R2 values
    ## (full/reduced) as well.
    mod00 <- matrix(rep(1, ncol(dat)), ncol = 1)
    n <- dim(dat)[2]
    m <- dim(dat)[1]
    df1 <- dim(mod)[2]
    df0 <- dim(mod0)[2]
    p <- rep(0, m)
    #Id <- diag(n)
    
    #resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
    resid <- dat - dat %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
    rss1 <- rowSums(resid * resid)
    rm(resid)
    
    if (df0 > 0)  {
        #resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
        resid0 <- dat - dat %*% mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
    } else {
        resid0 <- dat
    }
    rss0 <- rowSums(resid0 * resid0)
    rm(resid0)
    
    #resid00 <- dat %*% (Id - mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00))
    resid00 <- dat - dat %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
    rss00 <- rowSums(resid00 * resid00)
    rm(resid00)
    
    r2_full <- 1 - rss1/rss00
    r2_reduced <- 1 - rss0/rss00
    
    #delta <- (apply(dat, 1, mean) * 0.01)^2
    delta <- 0
    fstats <- ((rss0 - rss1)/(df1 - df0))/(delta + rss1/(n - df1))
    p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
    return(list(p = p, r2_full = r2_full, r2_reduced = r2_reduced))
}

batchEffectPvalue <- function(data, batch) {
    batch1 <- as.factor(batch)
    batch2 <- split(which(batch == batch1), batch1)
    meanbatch <- unlist(lapply(1:length(batch2), function(x) apply(data[, 
        batch2[[x]]], 1, mean)))
    varbatch <- unlist(lapply(1:length(batch2), function(x) apply(data[, 
        batch2[[x]]], 1, var)))
    skewbatch <- unlist(lapply(1:length(batch2), function(x) apply(data[, 
        batch2[[x]]], 1, skewness)))
    kurtbatch <- unlist(lapply(1:length(batch2), function(x) apply(data[, 
        batch2[[x]]], 1, kurtosis)))
    batcheffectmatrix <- rbind(meanbatch, varbatch, skewbatch, 
        kurtbatch)
    numgenes <- dim(data)[1]
    batch3 <- rep(1:length(batch2), each = numgenes)
    #genes <- rep(1:numgenes, length(batch2))
    #genes_mod <- model.matrix(~as.factor(genes))
    batch_mod <- model.matrix(~as.factor(batch3))
    #mod <- cbind(genes_mod, batch_mod[, -1])
    pdata <- data.frame(batch3)
    #mod00 <- model.matrix(~-1+0, data=pdata) 
    mod00 <- model.matrix(~1, data=pdata)
        # reduced model for batch adjusted data
    #batch_test <- delta_f.pvalue(batcheffectmatrix, mod, genes_mod)
    batch_test <- delta_f.pvalue(batcheffectmatrix, batch_mod, mod00)
    batch_ps <- batch_test$p
    return(batch_ps)
} 
