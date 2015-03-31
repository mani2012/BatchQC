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
biocLite(c('limma', 'preprocessCore', 'sva'))
```

Next, use [devtools](https://github.com/hadley/devtools) to install the latest
version of BatchQC from Github:
```r
require(devtools)
install_github("mani2012/BatchQC", build_vignettes=TRUE)
```

Install 'pandoc' package by following the instructions at the following URL:
http://johnmacfarlane.net/pandoc/installing.html

Install 'pander' package by following the instructions at the following URL:
http://rapporter.github.io/pander/

```r
require(devtools)
install_github('Rapporter/pander')
```

If all went well you should now be able to load BatchQC:
```r
require(BatchQC)
vignette('BatchQCIntro', package='BatchQC')
```

## Troubleshooting with Installation

If you are having issues with the installation, you may have to setup local directory, if you do not have permissions to install in the default location for R.
```r
export R_LIBS="/my_own_local_directory/R_libs"
```

And do something like the following
```r
install.packages("devtools", repos="http://cran.r-project.org", lib="/my_own_local_directory/R_libs")
```