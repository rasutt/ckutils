test_that("plot_exp_pop causes no errors", {
  # Try to plot expected population size
  plot_exp_pop(
    sim_yrs = 1:20, exp_N_t = 20*1.05^(1:20), ref_yr = 1, exp_N_ref = 20,
    srvy_yrs = 20
  )
})
