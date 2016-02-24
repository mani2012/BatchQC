BatchQC: Batch Effects Quality Control
======================================

The purpose of this package is to provide Quality Control of sequencing samples by deducing if there is batch effect and adjusts for it.

The package includes:

1. principal component analysis and plots to check batch effects
2. heatmap plot of gene expressions
3. statistical tests to assess batch effects
4. functions to generate simulated RNA-Seq data

`batchQC` is the pipeline function that generates the BatchQC report. It combines all the functions into one step.

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) along with a
few dependencies that BatchQC uses:

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('MCMCpack', 'limma', 'preprocessCore', 'sva', 'devtools', 'corpcor', 'matrixStats', 'shiny', 'ggvis', 'd3heatmap', 'reshape2', 'nlme', 'lme4', 'moments'))
```
Install 'pandoc' package by following the instructions at the following URL:
http://johnmacfarlane.net/pandoc/installing.html

Rstudio provide pandoc binaries at the following location for (windows, linux and mac):
https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip 

Install 'pander' package by following the instructions at the following URL:
http://rapporter.github.io/pander/

```r
require(devtools)
install_github('Rapporter/pander')
```
Install the latest development version of sva:
```r
require(devtools)
install_github('jtleek/sva-devel')
```

Next, use [devtools](https://github.com/hadley/devtools) to install the latest
version of HTShape and BatchQC from Github:
```r
require(devtools)
install_github("kokrah/HTShape")
install_github("mani2012/BatchQC", build_vignettes=TRUE)
```

If all went well you should now be able to load BatchQC:
```r
require(BatchQC)
vignette('BatchQCIntro', package='BatchQC')
vignette('BatchQC_examples', package='BatchQC')
```

## Troubleshooting with Installation

If you are having issues with the installation, you may have to setup local directory, if you do not have permissions to install in the default location for R. You may also want to load a version of R 3.1 or higher.
```r
export R_LIBS="/my_own_local_directory/R_libs"
module load R/R-3.1.1
```

And do something like the following
```r
install.packages("devtools", repos="http://cran.r-project.org", lib="/my_own_local_directory/R_libs")
```