
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
population and a sampling-survey study of it.

### Simulation years and expected population size

Close-kin studies use the probabilities of close family relationships
between pairs of sampled animals (close-kin pairs). These relationships
are defined by the breeding events that occur prior to sampling, so
ckutils lets you specify the number of years to simulate the population
prior to the final survey year. This needs to be long enough that most
of the relevant breeding events occur during the simulation.

ckutils simulates populations in which the expected population size
grows by a constant factor called the population growth rate $\lambda$.
Specifying this rate along with the expected population size $N$ in a
reference year defines the expected population size in all years. It is
good to check these as very small population sizes can be unsustainable,
and very large ones slow to simulate.

``` r
# Load close-kin utilities package
library(ckutils)

# Survey years and simulation length
srvy_yrs = c(2010, 2015, 2020)
sim_len = 50

# Population growth rate, and population size in reference year
lambda = 1.05
ref_yr = 2010
exp_N_ref = 1500

# Simulation years
fnl_yr = tail(srvy_yrs, 1)
init_yr = fnl_yr - sim_len + 1
sim_yrs = init_yr:fnl_yr

# Expected population size over simulation years
exp_N_init = exp_N_ref * lambda^(init_yr - ref_yr)
exp_N_t = exp_N_init * lambda^(0:(sim_len - 1))

# Plot expected population size over simulation years
plot_exp_pop(sim_yrs, exp_N_t, ref_yr, exp_N_ref, srvy_yrs)
```

<img src="man/figures/README-example-1.png" width="100%" />

### Predicted numbers of close-kin pairs

ckutils simulates populations in which each animal survives from one
year to the next with the same constant probability $\phi$, and becomes
sexually mature at the same age $\alpha$. Each animal also has the same
constant probability $p$ of being sampled in any survey (given that it
is alive in that year).

These parameters, along with the population and survey parameters above,
allow for approximations of the expected numbers of various close-kin
pairs among sampled animals. A sufficient number of such pairs is
essential to a successful study.

``` r
# Annual survival rate and age of sexual maturity
phi = 0.9
alpha = 8

# Capture probability
p = 0.1

# Survey year indices
s_yr_inds = seq_along(srvy_yrs)

# Number of surveys
k = length(srvy_yrs)

# Per capita birth rate
rho = lambda - phi

# Find expected numbers of kin-pairs in population
find_exp_ns_kps(exp_N_t, s_yr_inds, phi, rho, lambda, alpha, srvy_yrs, k)
#> $wtn
#>       N.s.yrs      APs     POPs     SMPs     SFPs     FSPs      HSPs
#> [1,] 223.7219 24913.89 265.6698 462.7737 477.7093 10.34484  919.7933
#> [2,] 234.9080 27473.44 278.9533 485.9124 501.6192 10.34484  966.8419
#> [3,] 246.6534 30295.64 292.9010 510.2080 526.7247 10.34484 1016.2430
#> 
#> $btn
#>           APs       SPs     POPs     SMPs SMPs.kwn.age     SFPs     FSPs
#> [1,] 52554.09 132.10557 566.1997 879.5848     4.159592 863.9176 17.60714
#> [2,] 55181.79  78.00702 545.6511 646.3162     2.020720 624.3719 11.86953
#> [3,] 57940.88 138.71085 594.5097 923.5640     4.367572 907.1424 17.60080
#>          HSPs
#> [1,] 1708.288
#> [2,] 1246.949
#> [3,] 1795.505
```

### Simulate population and study

``` r
# Number of loci in genome
L = 100

# Initial minor allele frequency
imaf = 0.5

# Set random seed for testing
set.seed(1)

# Simulate one population and study
pop_study = sim_pop_study(
  phi, lambda, exp_N_init, sim_len, srvy_yrs, k,
  fnl_yr, p, L, imaf, clvng.p = 0, tmp.emgn = 0,
  alpha = alpha, clvng.ints = F
)

# Look at it
head(pop_study)
#>        ID mum dad C2010 C2015 C2020 Cvg2010 Cvg2015 Cvg2020
#> 455   455 159  71     0     1     0       1       1       1
#> 506   506 115 150     1     0     0       0       0       0
#> 536   536 116 192     1     0     0       0       0       0
#> 1012 1012 411  65     0     1     0       0       0       0
#> 1035 1035  16 369     0     1     1       1       1       1
#> 1052 1052 478 480     1     0     0       0       0       0
```

``` r
names(attributes(pop_study))
#>  [1] "names"       "row.names"   "class"       "avg.phi.obs" "beta"       
#>  [6] "N.t.vec"     "ns.caps"     "Ns"          "ns.clvng"    "alive"      
#> [11] "alv.s.yrs"   "f.age"       "mum"         "dad"         "ID"         
#> [16] "ind.gts"
```
