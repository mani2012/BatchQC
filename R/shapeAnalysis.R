#' Perform Mean and Variance batch variation analysis
#' 
#' @param data Given data
#' @param groups a character vector indicating sample group membership
#' @param plot Indicate whether to generate plot
#' @param groupCol group color
#' @export
batchQC_shapeVariation = function (data, groups, plot=FALSE, groupCol=NULL) {
  
  sfit = fitMeanVariance(data)
  Y = sfit
  
  mod1 = model.matrix(~as.factor(groups), data=data.frame(Y)) # model with batch
  mod0 = model.matrix(~1, data=data.frame(Y)) # reduced model 
  
  ### Standard F test  ###
  n <- dim(Y)[1]
  Id <- diag(n)
  df1 <- dim(mod1)[2]
  df0 <- dim(mod0)[2]
  
  resid <- t(Y) %*% (Id - mod1 %*% solve(t(mod1) %*% mod1) %*% t(mod1))  # residuals full model
  rss1 <- rowSums(resid * resid)  ## SSE full model
  
  resid0 <- t(Y) %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)) # residuals reduced model
  rss0 <- rowSums(resid0 * resid0)   ## SSE reduced model
  
  Fstat = ( (rss0-rss1)/(df1-df0) ) / (rss1 /(n-df1) ) 
  # Compute p-value
  ps = 1-pf(Fstat,df1-df0,n-df1)
  mpval = ps[1]
  vpval = ps[2]
  mpvaltext <- round(mpval, 4)
  if (mpvaltext==0)  {
    mpvaltext <- "< 0.0001"
  }
  vpvaltext <- round(vpval, 4)
  if (vpvaltext==0)  {
    vpvaltext <- "< 0.0001"
  }
  
  if (plot) {
  
    pal = colorRampPalette(c("red", "orange", "white", "steelblue3", "navy"))(n=99)
    
    main = paste0("\nBatch Variation Analysis \n Mean Variation P-value = ", mpvaltext, 
                  "\n Variance Variation P-value = ", vpvaltext)
    
    gplots::heatmap.2(t(Y), trace="none", Rowv=F, Colv=F, dendrogram="none", col=pal,
                      ColSideColors=groupCol, density.info="none", scale="row",
                      cexRow=1.15,
                      colsep=cumsum(table(groups)), main=main)
    
    legend("bottomleft", legend=unique(groups), pch=19, col=unique(groupCol), 
           title="Batch")
  }
  ps
}

fitMeanVariance <- function (x) {
  # Compute Mean and Variance
  smean <- apply(x, 2, mean)
  svariance <- apply(x, 2, var)
  sfit <- cbind(smean, svariance)
  colnames(sfit) = c("Mean", "Variance")
  rownames(sfit) = colnames(x)
  return(sfit)
}
