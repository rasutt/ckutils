
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
good to check them for the entire simulation, as very small population
sizes can be unsustainable, and very large ones slow to simulate.

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

# Sampling probability
p = 0.1

# Survey year indices
s_yr_inds = seq_along(srvy_yrs)

# Number of surveys
k = length(srvy_yrs)

# Per capita birth rate
rho = lambda - phi

# Find expected numbers of kin-pairs in population
exp_ns_kps = find_exp_ns_kps(
  exp_N_t, s_yr_inds, phi, rho, lambda, alpha, srvy_yrs, k
)

# Combine expected numbers of kin-pairs within and between surveys in one matrix
exp_ns_kps_cmbd = cmbn_exp_ns_kps(exp_ns_kps, k)

# Find expected numbers of kin-pairs between sampled animals
exp_ns_kps_smpd = find_exp_ns_kps_smpd(exp_ns_kps_cmbd, k, p)
  
# Make it a data frame with row and column-names
exp_ns_kps_smpd_df = make_exp_ns_kps_smpd_df(srvy_yrs, exp_ns_kps_smpd)

# Display it nicely
knitr::kable(
  exp_ns_kps_smpd_df, 
  caption = "Predicted numbers of kin-pairs between sampled animals",
)
```

|                        |  2010 |  2015 |  2020 | 2010-2015 | 2010-2020 | 2015-2020 |
|:-----------------------|------:|------:|------:|----------:|----------:|----------:|
| Number sampled         |  22.4 |  23.5 |  24.7 |        NA |        NA |        NA |
| All pairs              | 239.1 | 264.2 | 291.9 |     525.5 |     551.8 |     579.4 |
| Self-pairs             |    NA |    NA |    NA |       1.3 |       0.8 |       1.4 |
| Parent-offspring pairs |   2.5 |   2.7 |   2.8 |       5.7 |       5.5 |       5.9 |
| Same-mother pairs      |   4.4 |   4.7 |   4.9 |       8.8 |       6.5 |       9.2 |
| Same-father pairs      |   4.6 |   4.8 |   5.1 |       8.6 |       6.2 |       9.1 |
| Full-sibling pairs     |   0.1 |   0.1 |   0.1 |       0.2 |       0.1 |       0.2 |
| Half-sibling pairs     |   8.8 |   9.3 |   9.8 |      17.1 |      12.5 |      18.0 |

Predicted numbers of kin-pairs between sampled animals

### Simulate population and study

ckutils simulates genetic inheritance of biallelic single nucleotide
polymorphisms (SNPs) from parents to offspring. The number of SNP loci
simulated $L$, and the initial minor allele frequency (IMAF) can be
specified.

The larger is $L$, the more informative the genotypes simulated, and the
better close-kin pairs can be distinguished from them, but the slower
the simulation and analysis is, and the more costly such genotypes are
to produce for real-world samples.

The closer the IMAF is to 0.5 the more informative each SNP locus is,
but in real-world genotypes MAFs vary. As the genetic inheritance
simulated in ckutils is probabilistic the MAFs at each locus will vary
independently, and at a rate that decreases with population size.

There are also options to specify biological scenarios such as
additional capture probability for females when they have a new
offspring, temporary emigration of males outside of sampling areas, and
selection of females to have offspring in order of time since last
having offspring.

``` r
# Number of loci in genome
L = 100

# Initial minor allele frequency
imaf = 0.5

# Set random seed for testing
set.seed(1)

# Simulate one population and study
pop_study = sim_pop_study(
  phi, lambda, exp_N_init, sim_len, srvy_yrs, k, fnl_yr, p, L, imaf, 
  clvng.p = 0, tmp.emgn = 0, alpha = alpha, clvng.ints = F
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
