% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_exp_pop.R
\name{plot_exp_pop}
\alias{plot_exp_pop}
\title{Plot expected population size over time}
\usage{
plot_exp_pop(sim_yrs, exp_N_t, ref_yr, exp_N_ref, srvy_yrs)
}
\arguments{
\item{sim_yrs}{Years over which population will be simulated}

\item{exp_N_t}{Expected population size in each year}

\item{ref_yr}{reference year at which expected population size is defined}

\item{exp_N_ref}{Expected population size in reference year}

\item{srvy_yrs}{Survey years}
}
\description{
Plot expected population size over time
}
\examples{
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

}
