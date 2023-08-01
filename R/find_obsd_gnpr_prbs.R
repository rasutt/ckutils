# Functions to find observed genopair probabilities given various kinships

#' Find sample genotypes, extracted from matrix of individual genotypes,
#' n_samples x n_loci, rows ordered by survey-year then individual ID
#'
#' @param smpl_hsts Sample histories. Binary matrix with rows for animals and
#'   columns for surveys.
#' @param indl_gts Individual genotypes (SNPs). Binary array with rows for
#'   maternal and paternal alleles, columns for genotype loci, and another
#'   dimension for animals.
#'
#' @return Binary array as for indl_gts but with the third dimension covering
#'   all samples, including repeated samples of the same animals.
#' @export
#'
#' @examples
find_smpl_gts <- function(smpl_hsts, indl_gts) {
  # Sample-individual indices, row numbers in sample history matrix,
  # representing the individual that each sample came from, ordered by
  # survey-year then individual ID
  smpl_indl_inds = row(smpl_hsts)[as.logical(smpl_hsts)]

  indl_gts[, , smpl_indl_inds]
}

#' Find log-probabilities of observed genopairs given that the individuals have
#' a particular kinship.  Computes over batches of loci to limit memory usage to
#' ~1Gb.
#'
#' @param pgps Possible genopair probabilities given various kinships,
#'   n_genotypes x n_genotypes x n_loci x n_kinships
#' @param ind_gts Individual genotypes. Genotypes of all individuals sampled
#' @param siips Sample-individual index pairs. The index of the individual that
#'   each sample comes from, for each sample in all pairs of samples.
#' @param L Number of genotype loci
#' @param sk Boolean indicating that only a single kinship is considered
#'
#' @return
#' @export
#'
#' @examples
find_oglps = function(pgps, ind_gts, siips, L, sk = F) {
  # Sample-individual indices of first and second samples in each pair, n_pairs
  siis_1 = siips[, 1]
  siis_2 = siips[, 2]

  # Number of pairs of samples
  n_pairs = length(siis_1)

  # Number of kinships for which genopair log-probabilities are being found
  n_knshps = if (sk) 1 else dim(pgps)[4]

  # Possible genopair log-probabilities at each locus given kinships,
  # n_genotypes x n_genotypes x n_loci x n_kinships
  pglps = log(pgps)

  # Actual genopair log-probabilities given kinships, n_pairs x n_kinships
  oglps = matrix(
    0, nrow = n_pairs, ncol = n_knshps,
    dimnames = list(
      pair = paste0("P", 1:n_pairs), kinship = dimnames(pgps)[[4]]
    )
  )

  # Genotype indices at each locus, n_individuals x n_loci, for genotypes
  # ordered 00, 01, 11, computed from binary alleles at each locus, 2 x n_loci x
  # n_individuals
  gt_inds = t(colSums(ind_gts)) + 1

  # Set size of batches of loci to keep memory usage < 1GB, n_pairs * btch_sz *
  # 8 bytes < 1Gb = 1e9 bytes ~ 2^30 => btch_sz ~< 2^27 / n_pairs
  btch_sz = 2^(26 - ceiling(log(n_pairs, 2)))

  # Find numbers of batches
  n_btchs = ceiling(L / btch_sz)

  # Set batch counter to zero
  btch_cnt = 0

  # Set timer
  s_time = proc.time()[3]

  # Loop over batches of loci
  for(btch_ind in 1:n_btchs) {
    # Increment batch counter
    btch_cnt = btch_cnt + 1

    # Find indices for current batch of loci
    loci_inds = ((btch_ind - 1) * btch_sz + 1):min(btch_ind * btch_sz, L)
    n_loci_btch = length(loci_inds)

    # Genotype indices for current batch of loci, n_individuals x n_loci_batch
    gt_inds_btch = gt_inds[, loci_inds]

    # Sample genotype indices for current batch of loci for first and second
    # samples in each pair, n_pairs x n_loci_batch
    sgisb_1 = gt_inds_btch[siis_1, ]
    sgisb_2 = gt_inds_btch[siis_2, ]

    # Sample genotype and locus-indices for current batch of loci for first and
    # second samples in each pair, (n_pairs x n_loci_batch) x 3
    sglisb = cbind(
      as.vector(sgisb_1),
      as.vector(sgisb_2),
      rep(loci_inds, each = n_pairs)
    )

    # Look up genopair log-probabilities for this batch of loci and add to
    # totals. Possible values are stored in n_gts x n_gts x L x kinships
    # array. Using apply to pick out one kinship at a time to limit memory
    # usage. No significant slow down as just a few kinships.
    if (sk) {
      oglps = oglps + rowSums(matrix(pglps[sglisb], nrow = n_pairs))
    } else {
      oglps = oglps +
        apply(pglps, 4, function(pglps_knshp) {
          rowSums(matrix(pglps_knshp[sglisb], nrow = n_pairs))
        })
    }
  }

  # Return genopair log-probabilities
  oglps
}

#' Find half-sibling versus unrelated pair PLODs
#'
#' @param oglps Observed genopair log-probabilities
#' @inheritParams find_oglps
#'
#' @return Vector with PLOD for each pair
#' @export
#'
#' @examples
find_hsp_up_plods <- function(oglps, L) {
  (oglps[, "HSP"] - oglps[, "UP"]) / L
}

#' Find expected values of half-sibling versus unrelated pair PLODs given
#' various kinships
#'
#' @inheritParams find_oglps
#'
#' @return Named vector
#' @export
#'
#' @examples
find_exp_plods <- function(pgps, L) {
  # Possible genopair log-probability ratios (given HSP vs UP) at each locus,
  # 3 x 3 x n_loci
  pglprs = log(pgps[, , , 2] / pgps[, , , 1])

  # Expected plods for kinship basis, unrelated, half-sibling, parent-offspring,
  # and self-pairs
  exp_plods_bss = colSums(pgps * rep(pglprs, 4), dims = 3) / L

  # Expected plods for extended kinships, first-cousin and avuncular pairs
  exp_plods_ext = (c(7, 3) * exp_plods_bss[1] + exp_plods_bss[3]) / c(8, 4)

  # Combine in order from furthest to closest kinship, add names, and return
  exp_plods = c(exp_plods_bss[1], exp_plods_ext, exp_plods_bss[2:4])
  names(exp_plods) = c(
    "Unrelated", "First cousin", "Avuncular", "Half-sibling",
    "Parent-offspring", "Self"
  )
  exp_plods
}

#' Plot half-sibling versus unrelated pair PLODs on log-scale
#'
#' @param hsp_up_plods Half-sibling versus unrelated pair PLODs
#' @param exp_plods Expected values of half-sibling versus unrelated pair PLODs
#'   for various kinships
#' @param show_legend Boolean indicating whether to show legend on plot
#'
#' @return None
#' @export
#'
#' @examples
plot_plods <- function(hsp_up_plods, exp_plods, show_legend = T) {
  # Get counts and take logs
  hst_data = graphics::hist(hsp_up_plods, breaks = 200, plot = F)
  hst_data$counts = log(hst_data$counts + 1)

  # Plot them
  plot(
    hst_data, main = "HSP vs UP PLODs for all samples",
    xlab = "PLOD", ylab = "Log (frequency + 1)"
  )

  # Add expected values
  graphics::abline(v = exp_plods, col = 2:7, lwd = 2)

  if (show_legend) {
    # Add legend
    graphics::legend(
      "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 0.5, lwd = 2,
      legend = c(
        "Expected value given kinship", "Unrelated", "First cousin",
        "Avuncular", "Half-sibling", "Parent-offspring", "Self"
      )
    )
  }
}

# Function to find genopair probabilities by exponentiating log-probabilities,
# checking for underflow, and trying to adjust if necessary
FindGPs = function(glps) {
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  gps = exp(glps)
  colnames(gps) = colnames(glps)
  all_undrflw = rowSums(gps) == 0

  # If there is underflow adjust log-probabilities by factor giving equal
  # weight to smallest and largest values to avoid both under and overflow
  if (any(all_undrflw)) {
    cat("Proportion of pairs for which probabilities given all kinships
        underflow to zero:", mean(all_undrflw), "\n")

    # Want smallest maximum kinship probability and largest probability to be
    # equally far from one
    adj = mean(c(min(apply(glps, 1, max)), max(glps)))
    glps_adj = glps - adj
    gps_adj = exp(glps_adj)
    colnames(gps_adj) = colnames(glps)

    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(glps_adj))
    print(summary(gps_adj))
    cat("Proportion of pairs for which adjusted probabilities given all
        kinships underflow to zero:", mean(rowSums(gps_adj) == 0), "\n")

    return(gps_adj)
  }
  else {
    return(gps)
  }
}

# # Function to find genopair probabilities from genotypes and sample-individual
# # index pairs
# FindGPsMdl <- function(pop.cap.hist, L, knshp.st, siips) {
#   # Get individual genotypes
#   gts = attributes(pop.cap.hist)$ind.gts
#
#   # Allele frequencies, 2 x n_loci matrices, representing relative
#   # frequencies of 0 and 1-coded SNP alleles at each locus
#   ale.frqs = FindAleFrqs(gts)
#
#   # Possible genopair probabilities given kinship set, n_possible_genotypes
#   # x n_possible_genotypes x n_loci x n_kinships, representing probabilities
#   # of each possible pair of genotypes at each locus given each kinship
#   # considered
#   pgps = FindPssGPPsKPs(ale.frqs, L, knshp.st)
#
#   # Genopair log-probabilities over all loci given each kinship, for each
#   # pair to include in likelihood
#   glps = FindGLPs(pgps, gts, siips, L)
#
#   # Exponentiate genopair log-probabilities given kinship set, checking for
#   # underflow and trying to adjust if necessary
#   FindGPs(glps)
# }
