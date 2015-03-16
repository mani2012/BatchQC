require(sva)
require(limma)
#' Checks for presence of batch effect and reports whether the batch needs to be adjusted
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @return pca Principal Components Analysis object of the data
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchQC_analyze(data.matrix, batch)
batchQC_analyze <- function(data.matrix, batch, mod=NULL)  {
  batchqc_heatmap(data.matrix, batch=batch, mod=mod)
  pca <- batchqc_pca(data.matrix, batch=batch, mod=mod)
  batchtest(pca, batch=batch, mod=mod)
  return(pca)
}

#' Checks for presence of batch effect and creates a html report 
#' with information including whether the batch needs to be adjusted
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param report_file Output report file name 
#' @param report_dir Output report directory path 
#' @return pca Principal Components Analysis object of the data
#' @export
#' @examples
#' nbatch <- 10
#' nperbatch <- 10
#' batch <- rep(1:nbatch, each=nperbatch)
#' batchQC(data.matrix, batch)
batchQC <- function(data.matrix, batch, mod=NULL, 
                           report_file="batchqc_report.html",
                           report_dir=".")  {
  rmdfile <- system.file("reports/batchqc_report.Rmd", package = "BatchQC")
  #rmarkdown::draft("batchqc_report.Rmd", template = "batchqc", package = "BatchQC")
  outputfile <- rmarkdown::render(rmdfile, output_file=report_file, output_dir=report_dir)
  return(outputfile)
}
