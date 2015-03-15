BatchQC: Batch Effects Quality Control
======================================

The purpose of this package is to provide Quality Control of sequencing samples by deducing if there is batch effect and adjusts for it.

The package includes:

1. principal component analysis and plots to check batch effects
2. heatmap plot of gene expressions
3. statistical tests to assess batch effects
4. functions to generate simulated RNA-Seq data

`batchQC` is the pipeline function. It combines all the functions into one step.
`batchQC_report` is the function that generates the BatchQC report.

## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) along with a
few dependencies that BatchQC uses:

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('limma', 'preprocessCore', 'sva'))
```

Next, use [devtools](https://github.com/hadley/devtools) to install the latest
version of cbcbSEQ from Github:
```r
require(devtools)
install_github("BatchQC", user="mani2012")
```

If all went well you should now be able to load BatchQC:
```r
require(BatchQC)
vignette('BatchQCIntro', package='BatchQC')
```

