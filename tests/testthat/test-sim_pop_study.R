test_that("SimPopStud() correct first row", {
  # Set random seed for testing
  set.seed(1)

  # Need a global variable setting births to be stochastic
  stch.bths = T

  # Simulate one population and study
  pop_study = sim_pop_study(
    phi = 0.9, lambda = 1.05, N.init = 20, hist.len = 20, srvy.yrs = 20, k = 1,
    f.year = 20, p = 0.5, L = 10, imaf = 0.5, clvng.p = 0, tmp.emgn = 0,
    alpha = 5, clvng.ints = F
  )

  expected = c(10, NA, NA, 1, 0)
  names(expected) = c('ID', 'mum', 'dad', 'C20', 'Cvg20')

  expect_equal(unlist(pop_study[1, ]), expected)
})
