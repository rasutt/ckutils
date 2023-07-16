test_that("find_exp_ns_kps works", {
  expect_snapshot(find_exp_ns_kps(
    exp.N.t = 20*1.05^(1:20), s.yr.inds = 1:2, phi = 0.9, rho = 0.15,
    lambda = 1.05, alpha = 5, srvy.yrs = c(18, 20), k = 2
  ))
})
