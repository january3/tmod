
<!-- README.md is generated from README.Rmd. Please edit that file -->
tmod
====

<!-- badges: start -->
<!-- badges: end -->
Tmod is a suite of gene set enrichment algorithms, visualizations and utilities which comes bundled with a few libraries of gene sets ("modules"). Following features distinguish tmod from other packages:

-   "panel plot" visualizations which allow to compare results of gene set enrichments;
-   several enrichment algorithms are implemented in tmod, in particular the very efficient, versatile and [reproducible](https://academic.oup.com/bioinformatics/article/35/24/5146/5511403) algorithm called "CERNO";
-   it includes a library of gene sets derived from clustering of gene expression data from human blood, which is especially useful in functional analysis in infection and immune responses.

An online version of tmod for demonstration purposes is available at [tmod.online](http://tmod.online).

Installation
------------

You can install the released version of tmod from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("tmod")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("january3/tmod")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(tmod)
#> For tmod user guide, type `vignette("tmod")`
## basic example code
```
