
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tmod <a href='https://january3.github.com/tmod/'><img src='man/figures/logo.png' align="right" height="138.5" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/january3/tmod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/january3/tmod/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Tmod is a suite of gene set enrichment algorithms, visualizations and
utilities which comes bundled with a few libraries of gene sets
(“modules”). Following features distinguish tmod from other packages:

-   “panel plot” visualizations which allow to compare results of gene
    set enrichments;
-   several enrichment algorithms are implemented in tmod, in particular
    the very efficient, versatile and
    [reproducible](https://academic.oup.com/bioinformatics/article/35/24/5146/5511403)
    algorithm called “CERNO”;
-   it includes a library of gene sets derived from clustering of gene
    expression data from human blood, which is especially useful in
    functional analysis in infection and immune responses.

## Installation

You can install the released version of tmod from
[CRAN](https://CRAN.R-project.org/) with:

``` r
install.packages("tmod")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("january3/tmod")
```

## Tmod manual

-   Access the documentation using `vignette("tmod")`
-   Full user manual is found [here](https://january3.github.io/tmod/)

## Example usage

``` r
library(tmod)
#> For tmod user guide, type `vignette("tmod")`
data(EgambiaResults)
tt <- EgambiaResults

## gene set enrichment analysis
res <- tmodCERNOtest(tt$GENE_SYMBOL)
head(res)
#>                ID                               Title    cerno  N1       AUC
#> LI.M37.0 LI.M37.0 immune activation - generic cluster 426.3578 100 0.7462103
#> DC.M4.2   DC.M4.2                        Inflammation 151.1520  20 0.9503953
#> DC.M3.4   DC.M3.4                          Interferon 129.4727  17 0.8315780
#> DC.M1.2   DC.M1.2                          Interferon 112.7056  17 0.9004196
#> DC.M7.29 DC.M7.29                        Undetermined 118.6759  20 0.8087599
#> LI.M11.0 LI.M11.0          enriched in monocytes (II) 113.8086  20 0.7766542
#>               cES      P.Value    adj.P.Val
#> LI.M37.0 2.131789 1.824844e-18 1.105856e-15
#> DC.M4.2  3.778799 8.040039e-15 2.436132e-12
#> DC.M3.4  3.808019 4.609405e-13 9.310998e-11
#> DC.M1.2  3.314869 2.298170e-10 3.481728e-08
#> DC.M7.29 2.966897 1.002268e-09 1.214749e-07
#> LI.M11.0 2.845216 5.255069e-09 5.307620e-07
```
