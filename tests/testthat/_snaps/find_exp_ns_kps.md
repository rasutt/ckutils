# find_exp_ns_kps works

    Code
      find_exp_ns_kps(exp.N.t = 20 * 1.05^(1:20), s.yr.inds = 1:2, phi = 0.9, rho = 0.15,
      lambda = 1.05, alpha = 5, srvy.yrs = c(18, 20), k = 2)
    Output
      $wtn
           N.s.yrs      APs     POPs     SMPs     SFPs     FSPs     HSPs
      [1,]   21.00 210.0000 24.93750 27.35514 28.00031 4.102449 47.15056
      [2,]   22.05 232.0763 26.18438 28.72290 29.41366 4.102449 49.93166
      
      $btn
              APs   SPs     POPs     SMPs SMPs.kwn.age     SFPs     FSPs     HSPs
      [1,] 463.05 17.01 52.09875 56.82054    0.3904464 56.61519 8.093546 97.24864
      

