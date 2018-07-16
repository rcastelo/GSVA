# GSVA: gene set variation analysis for microarray and RNA-seq data

[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/GSVA.svg)](http://bioconductor.org/packages/release/bioc/html/GSVA.html "How long has been GSVA in a release of Bioconductor")
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/GSVA.svg)](http://bioconductor.org/packages/stats/bioc/GSVA/ "Percentile (top 5/20/50% or 'available') of downloads over the last 6 full months")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/GSVA.svg)](http://bioconductor.org/packages/devel/bioc/html/GSVA.html#svn_source "Average SVN commits (to the devel branch) per month over the last 6 months")
[![Support posts](http://bioconductor.org/shields/posts/GSVA.svg)](https://support.bioconductor.org/t/GSVA/ "Bioconductor support site activity on GSVA, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")
<img align="right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/GSVA/GSVA.png" height="200"/>

**Current build status**
- `release` [![Bioconductor Availability](http://bioconductor.org/shields/availability/release/GSVA.svg)](http://bioconductor.org/packages/release/bioc/html/GSVA.html#archives "Whether GSVA release is available on all platforms") 
[![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/GSVA.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/GSVA/ "Bioconductor release build")
- `development` [![Bioconductor Availability](http://bioconductor.org/shields/availability/devel/GSVA.svg)](http://bioconductor.org/packages/devel/bioc/html/GSVA.html#archives "Whether GSVA devel is available on all platforms") 
[![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/GSVA.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/GSVA/ "Bioconductor devel build")

The `GSVA` package allows one to perform a change in coordinate systems of molecular measurements, transforming the data from a gene by sample matrix to a gene-set by sample matrix, thereby allowing the evaluation of pathway enrichment for each sample. This new matrix of GSVA enrichment scores facilitates applying standard analytical methods like functional enrichment, survival analysis, clustering, CNV-pathway analysis or cross-tissue pathway analysis, in a pathway-centric manner. For citing `GSVA` as a software package, please use the following reference:

  H&auml;nzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC _Bioinformatics_, 14:7, 2013.

## Installation

This is the __development__ version of the R/Bioconductor package GSVA. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [http://bioconductor.org/packages/GSVA](http://bioconductor.org/packages/GSVA) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the development version of R that you can find at [http://cran.r-project.org](http://cran.r-project.org) and then type the following instructions from the R shell:

```r
install.packages("BiocManager")
BiocManager::install("BiocInstaller", version="devel")
useDevel()
BiocManager::install("GSVA")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("rcastelo/GSVA")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **GSVA**
please use the [Bioconductor support site](http://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **GSVA**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/GSVA/issues](https://github.com/rcastelo/GSVA/issues)).
