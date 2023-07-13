
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ckutils

<!-- badges: start -->
<!-- badges: end -->

The goal of ckutils is to provide some utilities for close-kin study
design and implementation.

## Installation

You can install the development version of ckutils from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rasutt/ckutils")
```

## Example

This is a basic example which shows you how to simulate a wildlife
population and a sampling-survey study of it:

``` r
library(ckutils)

# Set random seed for testing
set.seed(1)

# Simulate one population and study
pop_stud = SimPopStud(
  phi = 0.9, lambda = 1.05, N.init = 20, hist.len = 20, srvy.yrs = 20, k = 1,
  f.year = 20, p = 0.5, L = 10, imaf = 0.5, clvng.p = 0, tmp.emgn = 0,
  alpha = 5, clvng.ints = F
)

# Look at it
head(pop_stud)
#>    ID mum dad C20 Cvg20
#> 10 10  NA  NA   1     0
#> 11 11  NA  NA   1     0
#> 20 20  NA  NA   1     0
#> 28 28  10  14   1     1
#> 32 32  21  18   1     0
#> 39 39   3  19   1     1
names(attributes(pop_stud))
#>  [1] "names"       "row.names"   "class"       "avg.phi.obs" "beta"       
#>  [6] "N.t.vec"     "ns.caps"     "Ns"          "ns.clvng"    "alive"      
#> [11] "alv.s.yrs"   "f.age"       "mum"         "dad"         "ID"         
#> [16] "ind.gts"
```
