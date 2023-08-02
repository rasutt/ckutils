#' Find simulation years
#'
#' @param fnl_yr Final year of simulation
#' @param sim_len Length of simulation
#'
#' @return Vector of years over which population will be simulated.
#' @export
#'
#' @examples
find_sim_yrs = function(fnl_yr, sim_len) {
  init_yr = fnl_yr - sim_len + 1
  init_yr:fnl_yr
}

find_exp_N_t <- function(exp_N_ref, lambda, init_yr, ref_yr, sim_len) {
  # Expected population size over simulation years
  exp_N_init = exp_N_ref * lambda^(init_yr - ref_yr)
  exp_N_t = exp_N_init * lambda^(0:(sim_len - 1))
}
