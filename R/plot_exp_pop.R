#' Plot expected population size over time
#'
#' @param sim_years Years over which population will be simulated
#' @param exp_N_t Expected population size in each year
#' @param base_yr Base year at which expected population size is determined
#' @param exp_N_base Expected population size in base year
#' @param srvy_yrs Survey years
#'
#' @export
#'
#' @examples
#' plot_exp_pop(
#'   sim_years = 1:20, exp_N_t = 20*1.05^(1:20), base_yr = 1, exp_N_base = 20,
#'   srvy_yrs = 20
#' )
#'
plot_exp_pop = function(sim_years, exp_N_t, base_yr, exp_N_base, srvy_yrs) {
  # Population size
  plot(
    sim_years, exp_N_t, ylim = c(0, max(exp_N_t)),
    col = 'red', lwd = 2, type = 'l',
    xlab = 'Year', ylab = 'Population size',
    main = "Expected population size over time"
  )

  # Base year and population
  graphics::abline(v = base_yr, h = exp_N_base, col = 2)

  # Surveys
  graphics::abline(v = srvy_yrs, lty = 2)

  # Add legend
  graphics::legend(
    "topleft", legend = c("Over time", "In base year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
}
