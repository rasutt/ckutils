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

# Function to find genopair log-probabilities given that the individuals have a
# particular kinship.  Computes over batches of loci to limit memory usage to
# ~1Gb.
FindGLPs = function(pgps, gts, siips, L, sk = F) {
  # Sample-individual indices of first and second samples in each pair, n_pairs
  siis.1 = siips[, 1]
  siis.2 = siips[, 2]

  # Number of pairs of samples
  n.pairs = length(siis.1)

  # Number of kinships for which genopair log-probabilities are being found
  n.knshps = if (sk) 1 else dim(pgps)[4]

  # Possible genopair log-probabilities at each locus given kinships,
  # n_genotypes x n_genotypes x n_loci x n_kinships
  pglps = log(pgps)

  # Actual genopair log-probabilities given kinships, n_pairs x n_kinships
  glps = matrix(
    0, nrow = n.pairs, ncol = n.knshps,
    dimnames = list(pair = NULL, Kinship = dimnames(pgps)[[4]])
  )

  # Genotype indices at each locus, n_individuals x n_loci, for genotypes
  # ordered 00, 01, 11, computed from binary alleles at each locus, 2 x n_loci x
  # n_individuals
  gt.inds = t(colSums(gts)) + 1

  # Set size of batches of loci to keep memory usage < 1GB, n.pairs * btch.sz *
  # 8 bytes < 1Gb = 1e9 bytes ~ 2^30 => btch.sz ~< 2^27 / n.pairs
  btch.sz = 2^(26 - ceiling(log(n.pairs, 2)))

  # Find numbers of batches
  n.btchs = ceiling(L / btch.sz)

  # Set batch counter to zero
  btch.cnt = 0

  # Set timer
  s.time = proc.time()[3]

  # Loop over batches of loci
  for(btch.ind in 1:n.btchs) {
    # Increment batch counter
    btch.cnt = btch.cnt + 1

    # Find indices for current batch of loci
    loci.inds = ((btch.ind - 1) * btch.sz + 1):min(btch.ind * btch.sz, L)
    n.loci.btch = length(loci.inds)

    # Genotype indices for current batch of loci, n_individuals x n_loci_batch
    gt.inds.btch = gt.inds[, loci.inds]

    # Sample genotype indices for current batch of loci for first and second
    # samples in each pair, n_pairs x n_loci_batch
    sgisb.1 = gt.inds.btch[siis.1, ]
    sgisb.2 = gt.inds.btch[siis.2, ]

    # Sample genotype and locus-indices for current batch of loci for first and
    # second samples in each pair, (n_pairs x n_loci_batch) x 3
    sglisb = cbind(
      as.vector(sgisb.1),
      as.vector(sgisb.2),
      rep(loci.inds, each = n.pairs)
    )

    # Look up genopair log-probabilities for this batch of loci and add to
    # totals. Possible values are stored in n_gts x n_gts x L x kinships
    # array. Using apply to pick out one kinship at a time to limit memory
    # usage. No significant slow down as just a few kinships.
    if (sk) {
      glps = glps + rowSums(matrix(pglps[sglisb], nrow = n.pairs))
    } else {
      glps = glps +
        apply(pglps, 4, function(pglps.knshp) {
          rowSums(matrix(pglps.knshp[sglisb], nrow = n.pairs))
        })
    }
  }

  # Return log genopair probabilities
  glps
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
    glps.adj = glps - adj
    gps.adj = exp(glps.adj)
    colnames(gps.adj) = colnames(glps)

    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(glps.adj))
    print(summary(gps.adj))
    cat("Proportion of pairs for which adjusted probabilities given all
        kinships underflow to zero:", mean(rowSums(gps.adj) == 0), "\n")

    return(gps.adj)
  }
  else {
    return(gps)
  }
}

# Function to find genopair probabilities from genotypes and sample-individual
# index pairs
FindGPsMdl <- function(pop.cap.hist, L, knshp.st, siips) {
  # Get individual genotypes
  gts = attributes(pop.cap.hist)$ind.gts

  # Allele frequencies, 2 x n_loci matrices, representing relative
  # frequencies of 0 and 1-coded SNP alleles at each locus
  ale.frqs = FindAleFrqs(gts)

  # Possible genopair probabilities given kinship set, n_possible_genotypes
  # x n_possible_genotypes x n_loci x n_kinships, representing probabilities
  # of each possible pair of genotypes at each locus given each kinship
  # considered
  pgps = FindPssGPPsKPs(ale.frqs, L, knshp.st)

  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  glps = FindGLPs(pgps, gts, siips, L)

  # Exponentiate genopair log-probabilities given kinship set, checking for
  # underflow and trying to adjust if necessary
  FindGPs(glps)
}
