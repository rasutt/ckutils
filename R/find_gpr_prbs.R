# Functions to find genopair probabilities given various kinships

#' Find allele frequencies
#'
#' @param ind_gts Genotypes in 2 x L x n_individuals arrays, representing two
#'   binary SNPs at each locus for each individual, excluding repeated samples
#'   of the same individual in different surveys
#'
#' @return 2 x L matrix representing the frequencies of 0 and 1-coded SNP
#'   alleles at each locus
#' @export
#'
#' @examples
find_ale_frqs <- function(ind_gts) {
  # Frequencies of 1-coded SNP alleles are means over both alleles at each
  # locus for each sample
  ale_frqs_1 = apply(ind_gts, 2, mean)

  # Combine with frequencies for 0-coded alleles and return
  ale_frqs = rbind(1 - ale_frqs_1, ale_frqs_1)
  dimnames(ale_frqs) =
    list(allele_code = paste0("A", 0:1), locus = colnames(ind_gts))
  ale_frqs
}

#' Find possible genotype probabilities
#'
#' @param ale_frqs 2 x L allele frequencies matrix for each allele of each
#'   possible genotype
#'
#' @return Number of possible genotypes x number of loci matrix
#' @export
#'
#' @examples
find_pss_gt_prbs = function(ale_frqs) {
  # Index and multiply the allele frequencies
  ales_1_inds = c(1, 1, 2)
  ales_2_inds = c(1, 2, 2)
  mat = ale_frqs[ales_1_inds, ] * ale_frqs[ales_2_inds, ]

  # Multiply by two possible cases for heterozygous genotypes
  mat[ales_1_inds != ales_2_inds, ] = mat[ales_1_inds != ales_2_inds, ] * 2
  dimnames(mat) = list(
    possible_genotype = c("GT00", "GT01", "GT11"), locus = colnames(ale_frqs)
  )
  mat
}


#' Find possible first genotype probabilities for genopairs, 3 x 3 x n_loci
#'
#' @param pss_gt_prbs Possible genotype probabilities
#' @param L Number of loci in each genotype
#'
#' @return
#' @export
#'
#' @examples
find_pss_frst_gt_prbs = function(pss_gt_prbs, L) {
  # 3 x L matrix indexed 3 times to fill the three columns.
  array(
    data = pss_gt_prbs[, rep(1:L, each = n_pss_gts)],
    dim = c(n_pss_gts, n_pss_gts, L),
    dimnames = pfgps_dimnames
  )
}

#' Find possible genopair probabilities for unrelated pairs
#'
#' @param pss_frst_gt_prbs Possible genotype probabilities for the first
#'   genotype in each genopair
#'
#' @return Array (Number of possible genotypes x number of possible genotypes x
#'   number of loci)
#' @export
#'
#' @examples
find_pss_gp_prbs_ups = function(pss_frst_gt_prbs) {
  # Conditional probabilities given that the pair are unrelated, the products of
  # the respective genotype probabilities, the second found by permuting the
  # array containing the first.
  pss_frst_gt_prbs * aperm(pss_frst_gt_prbs, c(2, 1, 3))
}

#' Find possible genopair probabilities for parent-offspring pairs
#'
#' @inheritParams find_pss_gt_prbs
#' @inheritParams find_pss_frst_gt_prbs
#' @inheritParams find_pss_gp_prbs_ups
#'
#' @return
#' @export
#'
#' @examples
find_pss_gp_prbs_pops = function(ale_frqs, L, pss_frst_gt_prbs) {
  # Conditional probabilities given that the pair are parent and offspring
  # (unordered). Products of the first genotype probabilities and the
  # conditional probabilities of the second genotypes given that the pair are
  # parent and offspring.
  pss_frst_gt_prbs *

    # Conditional probabilities of second genotype given that the pair are
    # parent and offspring (unordered).  0.5 for each allele in first genotype
    # being inherited, multiplied by the probability of the second genotype in
    # each case, as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin
    # Mark-Recapture. Filled into an L x 3 x 3 array which is permuted to the
    # standard dimensions.
    aperm(
      array(
        # Note: order data enters array is down columns, not across rows
        c(
          ale_frqs[1, ], 0.5 * ale_frqs[1, ], rep(0, L),
          ale_frqs[2, ], 0.5 * colSums(ale_frqs), ale_frqs[1, ],
          rep(0, L), 0.5 * ale_frqs[2, ], ale_frqs[2, ]
        ),
        c(L, n_pss_gts, n_pss_gts)
      ),
      c(2, 3, 1)
    )
}

#' Find possible genopair probabilities for self-pairs
#'
#' @inheritParams find_pss_gp_prbs_pops
#'
#' @return
#' @export
#'
#' @examples
find_pss_gp_prbs_sps = function(pss_frst_gt_prbs, L) {
  # Conditional probabilities given that the pair are the same individual
  # sampled twice. Genotype probabilities when the genotypes are the same, and
  # zero otherwise.
  pss_gp_prbs_sps = aperm(
    array(
      data = cbind(
        pss_frst_gt_prbs[1, 1, ], 0, 0,
        0, pss_frst_gt_prbs[2, 2, ], 0,
        0, 0, pss_frst_gt_prbs[3, 3, ]
      ),
      dim = c(L, n_pss_gts, n_pss_gts)
    ),
    c(2, 3, 1)
  )
  dimnames(pss_gp_prbs_sps) = pfgps_dimnames
  pss_gp_prbs_sps
}

#' Find possible genopair probabilities for half-sibling pairs
#'
#' @param pss_gp_prbs_ups
#' @param pss_gp_prbs_pops
#'
#' @return
#' @export
#'
#' @examples
find_pss_gp_prbs_hsps = function(pss_gp_prbs_ups, pss_gp_prbs_pops) {
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  (pss_gp_prbs_ups + pss_gp_prbs_pops) / 2
}

# Function to find possible genotype probabilities for all studies
FindPssGtPrbsAry = function(ale_frqs.ary) {
  # Indexing the 2 x L x n_sims allele frequencies array for each allele of
  # each possible genotype (globally defined for SNP genotypes), and multiplying
  # by 2 possible cases for heterozygous genotypes
  ary = ale_frqs.ary[ales.1.inds, , , drop = F] *
    ale_frqs.ary[ales.2.inds, , , drop = F]
  ary[ales.1.inds != ales.2.inds, , ] =
    ary[ales.1.inds != ales.2.inds, , , drop = F] * 2
  ary
}

# Function to find possible first genotype probabilities for genopairs for all
# studies, 3 x 3 x n_loci x n_sims
FindPssFrstGtPrbsAry = function(pss.gt.prbs.ary, L, n.sims) {
  # 3 x n_loci x n_sims array indexed 3 times to fill three columns.
  array(
    pss.gt.prbs.ary[, rep(1:L, each = n.pss.gts), ],
    c(n.pss.gts, n.pss.gts, L, n.sims)
  )
}

# Function to find possible genopair probabilities for unrelated pairs
FindPssGpPsUPsAry = function(pss.gt.1.prbs.ary) {
  # Products of the respective genotype probabilities, the second found by
  # permuting the array containing the first.
  pss.gt.1.prbs.ary * aperm(pss.gt.1.prbs.ary, c(2, 1, 3, 4))
}

# Function to find possible genopair probabilities for parent-offspring pairs
FindPssGpPsPOPsAry = function(pss.gt.1.prbs.ary, ale_frqs.ary, L, n.sims) {
  # Conditional probabilities given that the pair are parent and offspring
  # (unordered). Products of the first genotype probabilities and the
  # conditional probabilities of the second genotypes given that the pair are
  # parent and offspring.
  pss.gt.1.prbs.ary *

    # Conditional probabilities of second genotype given that the pair are
    # parent and offspring (unordered).  0.5 for each allele in first genotype
    # being inherited, multiplied by the probability of the second genotype in
    # each case, as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin
    # Mark-Recapture. Filled into an L x n_sims x 3 x 3 array which is permuted
    # to the standard dimensions.
    aperm(
      array(
        # Note: order data enters array is down columns, not across rows
        c(
          ale_frqs.ary[1, , ], 0.5 * ale_frqs.ary[1, , ], rep(0, L * n.sims),
          ale_frqs.ary[2, , ], rep(0.5, L * n.sims), ale_frqs.ary[1, , ],
          rep(0, L * n.sims), 0.5 * ale_frqs.ary[2, , ], ale_frqs.ary[2, , ]
        ),
        c(L, n.sims, n.pss.gts, n.pss.gts)
      ),
      c(3, 4, 1, 2)
    )
}

# Function to find possible genopair probabilities for self-pairs
FindPssGpPsSPsAry = function(pss.gt.prbs.ary, L, n.sims) {
  # Conditional probabilities given that the pair are the same individual
  # sampled twice. Genotype probabilities when the genotypes are the same, and
  # zero otherwise.
  aperm(
    array(
      cbind(
        as.vector(pss.gt.prbs.ary[1, , ]), 0, 0,
        0, as.vector(pss.gt.prbs.ary[2, , ]), 0,
        0, 0, as.vector(pss.gt.prbs.ary[3, , ])
      ),
      c(L, n.sims, n.pss.gts, n.pss.gts)
    ),
    c(3, 4, 1, 2)
  )
}

# Function to find possible genopair probabilities for half-sibling pairs
FindPssGpPsHSPs = function(pss.gp.prbs.UP, pss.gp.prbs.POP) {
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  (pss.gp.prbs.UP + pss.gp.prbs.POP) / 2
}

# Function to find possible genopair probabilities for half-sibling pairs
FindPssGpPsHSPsAry = function(pss.gp.prbs.UP.ary, pss.gp.prbs.POP.ary) {
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  (pss.gp.prbs.UP.ary + pss.gp.prbs.POP.ary) / 2
}

# Function to find possible genopair probabilities over multiple loci given
# multiple kinships
FindPssGPPsKPs = function(ale_frqs, L, knshp.st) {
  # Find possible genotype probabilities
  pss.gt.prbs = FindPssGtPrbs(ale_frqs)

  # Find possible first genotype probabilities for genopairs
  pss.gt.1.prbs = FindPssFrstGtPrbs(pss.gt.prbs, L)

  # Find possible genopair probabilities for unrelated pairs
  pss.gp.prbs.UP = FindPssGpPsUPs(pss.gt.1.prbs)

  # If parent-offspring pairs in kinship set selected
  if("Parent-offspring" %in% knshp.st) {
    # Find possible genopair probabilities for parent-offspring pairs
    pss.gp.prbs.POP = FindPssGpPsPOPs(ale_frqs, L, pss.gt.1.prbs)
  } else {
    pss.gp.prbs.POP = NULL
  }

  # If self-pairs in kinship set selected
  if("Self" %in% knshp.st) {
    # Function to find possible genopair probabilities for self-pairs
    pss.gp.prbs.SP = FindPssGpPsSPs(pss.gt.1.prbs, L)
  } else {
    pss.gp.prbs.SP = NULL
  }

  # If half-sibling pairs in kinship set selected
  if("Half-sibling" %in% knshp.st) {
    # Find possible genopair probabilities for half-sibling pairs
    pss.gp.prbs.HSP = FindPssGpPsHSPs(pss.gp.prbs.UP, pss.gp.prbs.POP)
  } else {
    pss.gp.prbs.HSP = NULL
  }

  # Combine genopair probabilities for selected kinships and return
  array(
    c(pss.gp.prbs.UP, pss.gp.prbs.HSP, pss.gp.prbs.POP, pss.gp.prbs.SP),
    dim = c(n.pss.gts, n.pss.gts, L, length(knshp.st) + 1),
    dimnames = list(
      gt.1 = pss.gt.lbls, gt.2 = pss.gt.lbls, Locus = paste0("L", 1:L),
      Kinship = c("UP", c("HSP", "POP", "SP")[rev(knshp.chcs) %in% knshp.st])
    )
  )
}
