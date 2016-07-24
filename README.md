BatchQC: Batch Effects Quality Control
======================================

Sequencing and microarray samples often are collected or processed in multiple 
batches or at different times. This often produces technical biases that can 
lead to incorrect results in the downstream analysis. BatchQC is a software tool
that streamlines batch preprocessing and evaluation by providing interactive 
diagnostics, visualizations, and statistical analyses to explore the extent to 
which batch variation impacts the data. BatchQC diagnostics help determine 
whether batch adjustment needs to be done, and how correction should be applied 
before proceeding with a downstream analysis. Moreover, BatchQCvinteractively 
applies multiple common batch effect approaches to the data, and the user can 
quickly see the benefits of each method. BatchQC is developed as a Shiny App. 
The output is organized into multiple tabs, and each tab features an important 
part of the batch effect analysis and visualization of the data. The BatchQC 
interface has the following analysis groups: Summary, Differential Expression, 
Median Correlations, Heatmaps, Circular Dendrogram, PCA Analysis, Shape, ComBat 
and SVA. 

The package includes:

1. Summary and Sample Diagnostics
2. Differential Expression Plots and Analysis using LIMMA
3. Principal Component Analysis and plots to check batch effects
4. Heatmap plot of gene expressions
5. Median Correlation Plot
6. Circular Dendrogram clustered and colored by batch and condition
7. Shape Analysis for the distribution curve based on HTShape package
8. Batch Adjustment using ComBat
9. Surrogate Variable Analysis using sva package
10. Function to generate simulated RNA-Seq data

`batchQC` is the pipeline function that generates the BatchQC report and 
launches the Shiny App in interactive mode. It combines all the functions into 
one step.

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) and simply
run the following to automatically install BatchQC and all the dependencies, 
except pandoc, which you have to manually install as follows.

```r
source("http://bioconductor.org/biocLite.R")
biocLite("BatchQC")
```
Install 'pandoc' package by following the instructions at the following URL:
http://pandoc.org/installing.html

Rstudio also provides pandoc binaries at the following location for Windows, 
Linux and Mac:
https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip 

If all went well you should now be able to load BatchQC:
```r
require(BatchQC)
browseVignettes("BatchQC")
vignette('BatchQCIntro', package='BatchQC')
vignette('BatchQC_examples', package='BatchQC')
vignette('BatchQC_usage_advanced', package='BatchQC')
```

If you want to install the latest development version of BatchQC from Github, 
Use [devtools](https://github.com/hadley/devtools) to install it as follows:
```r
require(devtools)
install_github("mani2012/BatchQC", build_vignettes=TRUE)
```
If you want to manually install the BatchQC dependencies, run the following:
```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('MCMCpack', 'limma', 'preprocessCore', 'sva', 'devtools', 'corpcor', 
'matrixStats', 'shiny', 'ggvis', 'd3heatmap', 'reshape2', 'moments', 
'rmarkdown', 'pander', 'gplots'))
```

## Troubleshooting with Installation

Please make sure you have installed pandoc by following the instructions from http://pandoc.org/installing.html. Otherwise, you may get an error such as the 
following:

```r
* creating vignettes ... ERROR
Warning in engine$weave(file, quiet = quiet, encoding = enc) :
  Pandoc (>= 1.12.3) and/or pandoc-citeproc is not available. Please install 
  both.
Error: processing vignette 'BatchQCIntro.Rmd' failed with diagnostics:
It seems you should call rmarkdown::render() instead of knitr::knit2html() 
because BatchQCIntro.Rmd appears to be an R Markdown v2 document.
Execution halted
Error: Command failed (1)
```
For generating pdf vignettes in Linux, you need to install texlive and lmodern 
as follows:

```r
sudo apt-get install texlive
sudo apt-get install lmodern
```

If you do not have permissions to install in the default location for R, you 
may have to setup local directory. You may also want to load a version of 
R 3.3.0 or higher.
```r
export R_LIBS="/my_own_local_directory/R_libs"
module load R/R-3.3.0
```

And do something like the following
```r
install.packages("devtools", repos="http://cran.r-project.org", 
lib="/my_own_local_directory/R_libs")
```

If you get an error "X11 font -adobe-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*, 
face 1 at size 6 could not be loaded" in your browser, when you view
some plots, this could be due to a missing fonts setup in the OS.
This issue has been discussed at the following link:
http://goo.gl/ukQXMI

The fonts installed can be found by typing the command ‘xlsfonts’ from bash and 
installing (or re-installing) xorg-x11-fonts for 75 and 100 dpi as mentioned
in the link above could do the trick to resolve this issue.
