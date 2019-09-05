# GSVA: gene set variation analysis for microarray and RNA-seq data

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/GSVA.svg)](https://bioconductor.org/packages/release/bioc/html/GSVA.html "How long has been GSVA in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/GSVA.svg)](https://bioconductor.org/packages/stats/bioc/GSVA/ "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/GSVA.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GSVA "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Support posts](https://bioconductor.org/shields/posts/GSVA.svg)](https://support.bioconductor.org/t/GSVA/ "Support site activity on GSVA, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")
<img align="right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/GSVA/GSVA.png" height="200"/>

**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/GSVA.svg)](https://bioconductor.org/packages/release/bioc/html/GSVA.html#archives "Whether GSVA release is available on all platforms") 
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/GSVA.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/GSVA/ "Bioconductor release build")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/GSVA.svg)](https://bioconductor.org/packages/release/bioc/html/GSVA.html#since "Number of recursive dependencies needed to install package")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/GSVA.svg)](https://bioconductor.org/packages/devel/bioc/html/GSVA.html#archives "Whether GSVA devel is available on all platforms") 
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/GSVA.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GSVA/ "Bioconductor devel build")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/GSVA.svg)](https://bioconductor.org/packages/devel/bioc/html/GSVA.html#since "Number of recursive dependencies needed to install package")

The `GSVA` package allows one to perform a change in coordinate systems of molecular measurements, transforming the data from a gene by sample matrix to a gene-set by sample matrix, thereby allowing the evaluation of pathway enrichment for each sample. This new matrix of GSVA enrichment scores facilitates applying standard analytical methods such as functional enrichment, survival analysis, clustering, CNV-pathway analysis or cross-tissue pathway analysis, in a pathway-centric manner. For citing `GSVA` as a software package, please use the following reference:

  H&auml;nzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC _Bioinformatics_, 14:7, 2013.

## Installation

This is the __development__ version of the R/Bioconductor package GSVA. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [https://bioconductor.org/packages/GSVA](https://bioconductor.org/packages/GSVA) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the development version of R that you can find at [https://cran.r-project.org](https://cran.r-project.org) and then type the following instructions from the R shell:

```r
install.packages("BiocManager")
BiocManager::install("GSVA", version = "devel")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("rcastelo/GSVA")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **GSVA**
please use the [Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **GSVA**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/GSVA/issues](https://github.com/rcastelo/GSVA/issues)).
