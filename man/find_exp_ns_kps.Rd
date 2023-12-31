% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_exp_ns_kps.R
\name{find_exp_ns_kps}
\alias{find_exp_ns_kps}
\title{Find expected numbers of kin-pairs in population}
\usage{
find_exp_ns_kps(exp.N.t, s.yr.inds, phi, rho, lambda, alpha, srvy.yrs, k)
}
\arguments{
\item{exp.N.t}{Expected population size over time}

\item{s.yr.inds}{Survey-year indices}

\item{phi}{Individual survival rate}

\item{rho}{Per capita birth rate}

\item{lambda}{Population growth rate}

\item{alpha}{Age of sexual maturity}

\item{srvy.yrs}{Survey-years}

\item{k}{Number of survey-years}
}
\value{
A list of two matrices, both with rows for kinship-types, one with
columns for survey-years and one for pairs of survey-years.
}
\description{
Find expected numbers of kin-pairs in population
}
\examples{
find_exp_ns_kps(
  exp.N.t = 20*1.05^(1:20), s.yr.inds = 1:2, phi = 0.9, rho = 0.15,
  lambda = 1.05, alpha = 5, srvy.yrs = c(18, 20), k = 2
)

}
