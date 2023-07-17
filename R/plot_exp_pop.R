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
#' plot_exp_pop(
#'   sim_yrs = 1:20, exp_N_t = 20*1.05^(1:20), ref_yr = 1, exp_N_ref = 20,
#'   srvy_yrs = 20
#' )
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
