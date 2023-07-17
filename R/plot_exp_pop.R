#' Plot expected population size over time
#'
#' @param sim_yrs Years over which population will be simulated
#' @param exp_N_t Expected population size in each year
#' @param ref_yr reference year at which expected population size is defined
#' @param exp_N_ref Expected population size in reference year
#' @param srvy_yrs Survey years
#'
#' @export
#'
#' @examples
#' # Survey years and simulation length
#' srvy_yrs = c(2010, 2015, 2020)
#' sim_len = 50
#'
#' # Population growth rate, and population size in reference year
#' lambda = 1.05
#' ref_yr = 2010
#' exp_N_ref = 1500
#'
#' # Simulation years
#' fnl_yr = tail(srvy_yrs, 1)
#' init_yr = fnl_yr - sim_len + 1
#' sim_yrs = init_yr:fnl_yr
#'
#' # Expected population size over simulation years
#' exp_N_init = exp_N_ref * lambda^(init_yr - ref_yr)
#' exp_N_t = exp_N_init * lambda^(0:(sim_len - 1))
#'
#' # Plot expected population size over simulation years
#' plot_exp_pop(sim_yrs, exp_N_t, ref_yr, exp_N_ref, srvy_yrs)
#'
plot_exp_pop = function(sim_yrs, exp_N_t, ref_yr, exp_N_ref, srvy_yrs) {
  # Population size
  plot(
    sim_yrs, exp_N_t, ylim = c(0, max(exp_N_t)),
    col = 'red', lwd = 2, type = 'l',
    xlab = 'Year', ylab = 'Population size',
    main = "Expected population size over simulation"
  )

  # ref year and population
  graphics::abline(v = ref_yr, h = exp_N_ref, col = 2)

  # Surveys
  graphics::abline(v = srvy_yrs, lty = 2)

  # Add legend
  graphics::legend(
    "topleft", legend = c("Over time", "In reference year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
}
