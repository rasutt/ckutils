#' Find expected numbers of kin-pairs in population
#'
#' @param exp.N.t Expected population size over time
#' @param s.yr.inds Survey-year indices
#' @param phi Individual survival rate
#' @param rho Per capita birth rate
#' @param lambda Population growth rate
#' @param alpha Age of sexual maturity
#' @param srvy.yrs Survey-years
#' @param k Number of survey-years
#'
#' @return A list of two matrices, both with rows for kinship-types, one with
#'   columns for survey-years and one for pairs of survey-years.
#' @export
#'
#' @examples
#' find_exp_ns_kps(
#'   exp.N.t = 20*1.05^(1:20), s.yr.inds = 1:2, phi = 0.9, rho = 0.15,
#'   lambda = 1.05, alpha = 5, srvy.yrs = c(18, 20), k = 2
#' )
#'
find_exp_ns_kps = function(
    exp.N.t, s.yr.inds, phi, rho, lambda, alpha, srvy.yrs, k
) {
  ## Intermediate results

  # Expected population sizes in survey years
  exp.N.s.yrs = exp.N.t[s.yr.inds]
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  # Probability not new-born (phi over lambda)
  p.o.l = phi / lambda
  # Reciprocal of probability not new-born (phi over lambda)
  l.o.p = lambda / phi
  # Reciprocal of probability that an animal is mature
  rcl.prb.mtr = l.o.p^alpha
  # Birth rate among mature females
  beta = 2 * (1 - p.o.l) * rcl.prb.mtr

  ### Expected numbers of kin-pairs

  ## Within surveys

  # All pairs in survey years
  exp.ns.APs.wtn = choose(exp.N.s.yrs, 2)

  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = exp.N.s.yrs * rho * (1 + phi) / lmb.m.ph.sq

  # Same-mother pairs within survey years
  exp.ns.SMPs.wtn = exp.N.s.yrs * beta * rho * phi^2 / lmb.m.ph.sq^2

  # Same-father pairs within survey years, split into same and different birth
  # years as used separately later
  exp.ns.SFPs.diff.b.yrs.wtn = phi * exp.ns.SMPs.wtn
  exp.ns.SFPs.same.b.yr.wtn = beta^2 * phi^(alpha + 1) / 4 *
    (exp.N.s.yrs / (lambda^(alpha - 1) * lmb.m.ph.sq) - 1 / (1 - phi^2))
  exp.ns.SFPs.wtn = exp.ns.SFPs.diff.b.yrs.wtn + exp.ns.SFPs.same.b.yr.wtn

  # Full-sibling pairs within survey years, constant over time but gets repeated
  # by cbind when returned
  exp.ns.FSPs.wtn = 2 * beta * rcl.prb.mtr * rho * phi^4 /
    (lambda * (lambda - phi^3) * (1 - phi^2))

  # Half-sibling pairs within survey years
  exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn

  ## Between surveys

  # Survey-pair indices
  s.pr.inds = utils::combn(k, 2)
  s.inds.1 = s.pr.inds[1, ]
  s.inds.2 = s.pr.inds[2, ]
  exp.N.s.yrs.1 = exp.N.s.yrs[s.inds.1]
  exp.N.s.yrs.2 = exp.N.s.yrs[s.inds.2]
  s.yrs.1 = srvy.yrs[s.inds.1]
  s.yrs.2 = srvy.yrs[s.inds.2]
  s.gaps = s.yrs.2 - s.yrs.1
  p.t.s.gs = phi^s.gaps

  # All pairs between pairs of surveys
  exp.ns.APs.btn = exp.N.s.yrs.1 * exp.N.s.yrs.2

  # Self-pairs between pairs of survey years
  exp.ns.SPs.btn = p.t.s.gs * exp.N.s.yrs.1

  # Parent-offspring pairs between survey years
  s.yrs.1.p.a = s.yrs.1 + alpha
  exp.ns.POPs.btn = 2 * exp.ns.POPs.wtn[s.inds.1] * p.t.s.gs +
    2 * exp.N.s.yrs.2 * (1 - p.o.l) * p.o.l^s.yrs.2 *
    ((l.o.p^(s.yrs.1 + 1) - l.o.p^(pmin(s.yrs.1.p.a, s.yrs.2) + 1)) /
       (1 - l.o.p) +
       ifelse(
         s.yrs.1.p.a < s.yrs.2,
         (s.yrs.2 - s.yrs.1.p.a) * l.o.p^s.yrs.1.p.a,
         0
       ))

  # Same-mother pairs between survey years
  exp.ns.SMPs.btn = 2 * exp.ns.SMPs.wtn[s.inds.1] * p.t.s.gs +
    s.gaps * exp.N.s.yrs.2 * beta * (1 - p.o.l) *
    lambda / lmb.m.ph.sq * p.o.l^s.gaps

  # Same-mother pairs between surveys, ages known (five and zero)
  exp.ns.SMPs.kwn.age.btn = exp.N.s.yrs.2 * (1 - p.o.l) * beta *
    phi^5 * p.o.l^(s.gaps + 5)

  # Same-father pairs between survey years
  exp.ns.SFPs.btn = phi * exp.ns.SMPs.btn +
    2 * p.t.s.gs * exp.ns.SFPs.same.b.yr.wtn[s.inds.1]

  # Full-sibling pairs between survey years, note the predicted number within
  # surveys is constant
  exp.ns.FSPs.btn = 2 * exp.ns.FSPs.wtn * p.t.s.gs +
    2 * beta * (1 - p.o.l) * rcl.prb.mtr * phi^(s.yrs.1 + s.yrs.2 + 1) *
    (p.o.l^(s.yrs.1 + 1) - p.o.l^(s.yrs.2 + 1)) /
    ((phi^3 / lambda)^s.yrs.1 * (1 - phi^3 / lambda) * (1 - p.o.l))

  # Half-sibling pairs between surveys
  exp.ns.HSPs.btn = exp.ns.SMPs.btn + exp.ns.SFPs.btn - 2 * exp.ns.FSPs.btn

  # Return as list
  list(
    wtn = cbind(
      N.s.yrs = exp.N.s.yrs, APs = exp.ns.APs.wtn, POPs = exp.ns.POPs.wtn,
      SMPs = exp.ns.SMPs.wtn, SFPs = exp.ns.SFPs.wtn, FSPs = exp.ns.FSPs.wtn,
      HSPs = exp.ns.HSPs.wtn
    ),
    btn = cbind(
      APs = exp.ns.APs.btn, SPs = exp.ns.SPs.btn, POPs = exp.ns.POPs.btn,
      SMPs = exp.ns.SMPs.btn, SMPs.kwn.age = exp.ns.SMPs.kwn.age.btn,
      SFPs = exp.ns.SFPs.btn, FSPs = exp.ns.FSPs.btn, HSPs = exp.ns.HSPs.btn
    )
  )
}

#' Combine expected numbers of kin-pairs within and between surveys in one
#' matrix
#'
#' @inheritParams find_exp_ns_kps
#' @param exp_ns_kps Expected numbers of kin-pairs, as output by
#'   find_exp_ns_kps. A list of two matrices, both with rows for survey-years
#'   and one for pairs of survey-years, one with columns for kinship-types
#'
#' @return A matrix with rows for survey years and pairs, and columns for
#'   kin-pair types
#' @export
#'
#' @examples
#'
#'
cmbn_exp_ns_kps = function(exp_ns_kps, k) {
  rbind(
    # Within surveys
    cbind(
      # Population sizes and total numbers of pairs
      exp_ns_kps$wtn[, 1:2],

      # Self-pairs don't apply
      rep(NA, k),

      # Other close-kin pairs
      exp_ns_kps$wtn[, -(1:2)]
    ),

    # Between survey-pairs
    cbind(
      # Population sizes don't apply
      rep(NA, choose(k, 2)),

      # Other close-kin pairs
      exp_ns_kps$btn[, -5]
    )
  )
}

#' Find expected numbers of kin-pairs between sampled animals
#'
#' @inheritParams find_exp_ns_kps
#' @param exp_ns_kps_cmbd Expected numbers of kin-pairs, as output by
#'   cmbn_exp_ns_kps. A matrix with rows for survey years and pairs, and columns
#'   for kin-pair types
#' @param p Sampling probability
#'
#' @return A matrix with rows for survey years and pairs, and columns for
#'   kin-pair types
#' @export
#'
#' @examples
#'
find_exp_ns_kps_smpd = function(exp_ns_kps_cmbd, k, p) {
  # Expected numbers of animals sampled
  exp_ns_smpd = exp_ns_kps_cmbd[1:k, "N.s.yrs"] * p

  # Predicted total numbers of pairs among sampled animals, within surveys, and
  # between survey-pairs
  pred_ns_APs_smpd = c(
    choose(exp_ns_smpd, 2),
    utils::combn(exp_ns_smpd, 2, function(ns) ns[1] * ns[2])
  )

  cbind(
    # Numbers sampled, not applicable to survey-pairs
    c(exp_ns_smpd, rep(NA, choose(k, 2))),

    # Total numbers of pairs
    pred_ns_APs_smpd,

    # Close-kin pairs. Probabilities from numbers in population multiplied by
    # total numbers of pairs.
    exp_ns_kps_cmbd[, -(1:2)] / exp_ns_kps_cmbd[, 2] * pred_ns_APs_smpd
  )
}

#' Make data frame for expected numbers of kin-pairs between sampled animals
#'
#' @inheritParams find_exp_ns_kps
#' @inheritParams find_exp_ns_kps_smpd
#'
#' @return A data frame with rows for survey years and pairs, and columns for
#'   kin-pair types
#' @export
#'
#' @examples
#'
make_exp_ns_kps_smpd_df <- function(
    exp_N_t, s_yr_inds, phi, rho, lambda, alpha, srvy_yrs, k, p
) {
  # Find expected numbers of kin-pairs in population
  exp_ns_kps = find_exp_ns_kps(
    exp_N_t, s_yr_inds, phi, rho, lambda, alpha, srvy_yrs, k
  )

  # Combine numbers within and between surveys in one matrix
  exp_ns_kps_cmbd = cmbn_exp_ns_kps(exp_ns_kps, k)

  # Find numbers among only sampled animals
  exp_ns_kps_smpd = find_exp_ns_kps_smpd(exp_ns_kps_cmbd, k, p)

  # Kin-pair types
  kpts = c(
    "Number sampled", "All pairs", "Self-pairs",
    "Parent-offspring pairs",
    "Same-mother pairs", "Same-father pairs", "Full-sibling pairs",
    "Half-sibling pairs"
  )

  # Survey year pairs
  s_yr_prs = apply(combn(srvy_yrs, 2), 2, paste, collapse = "-")

  # Round and make data frame with row and column-names
  exp_ns_kps_smpd_df = data.frame(
    round(t(exp_ns_kps_smpd), 1), row.names = kpts
  )
  names(exp_ns_kps_smpd_df) = c(srvy_yrs, s_yr_prs)

  # Return it
  exp_ns_kps_smpd_df
}
