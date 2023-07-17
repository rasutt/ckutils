test_that("plot_exp_pop causes no errors", {
  # Try to plot expected population size
  expect_no_error({
    # Population size in reference year
    ref_yr = 2010
    exp_N_ref = 1500

    # Population growth rate
    lambda = 1.05

    # Survey years
    srvy_yrs = c(2010, 2015, 2020)

    # Simulation length
    sim_len = 50

    # Simulation years
    fnl_yr = tail(srvy_yrs, 1)
    init_yr = fnl_yr - sim_len + 1
    sim_yrs = init_yr:fnl_yr

    # Expected population size over simulation years
    exp_N_init = exp_N_ref * lambda^(init_yr - ref_yr)
    exp_N_t = exp_N_init * lambda^(0:(sim_len - 1))

    # Plot expected population size over simulation years
    plot_exp_pop(sim_yrs, exp_N_t, ref_yr, exp_N_ref, srvy_yrs)
  })
})
