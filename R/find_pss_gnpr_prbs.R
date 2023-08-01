# Functions to find possible genopair probabilities given various kinships

# Set number of possible genotypes to allow for more later
n_pss_gts = 3

# Set genotype names
gt_names = c("GT00", "GT01", "GT11")

# Set genopair matrix row and column-names
gnpr_names = list(first_genotype = gt_names, second_genotype = gt_names)

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
    possible_genotype = gt_names, locus = colnames(ale_frqs)
  )
  mat
}

#' Round numeric matrix to two decimal places and convert to data frame for
#' display
#'
#' @param mat Numeric matrix
#'
#' @return Data frame rounded to two decimal places
#' @export
#'
#' @examples
dfr2 = function(mat) {
  data.frame(round(mat, 2))
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
    dimnames =  c(gnpr_names, list(locus = colnames(pss_gt_prbs)))
  )
}

#' Find possible genopair probabilities for unrelated pairs
#'
#' @param pss_frst_gt_prbs Possible genotype probabilities for the first
#'   genotype in each genopair
#'
#' @return Array with dimensions (number of possible genotypes x number of
#'   possible genotypes x number of loci)
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
#' @return Array with dimensions (number of possible genotypes x number of
#'   possible genotypes x number of loci)
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
#' @return Array with dimensions (number of possible genotypes x number of
#'   possible genotypes x number of loci)
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
  dimnames(pss_gp_prbs_sps) = dimnames(pss_frst_gt_prbs)
  pss_gp_prbs_sps
}

#' Find possible genopair probabilities for half-sibling pairs
#'
#' @param pss_gp_prbs_ups Possible genopair probabilities for unrelated pairs
#' @param pss_gp_prbs_pops Possible genopair probabilities for parent-offspring
#'   pairs
#'
#' @return Array with dimensions (number of possible genotypes x number of
#'   possible genotypes x number of loci)
#' @export
#'
#' @examples
find_pss_gp_prbs_hsps = function(pss_gp_prbs_ups, pss_gp_prbs_pops) {
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  (pss_gp_prbs_ups + pss_gp_prbs_pops) / 2
}

#' Find possible genopair probabilities over multiple loci given multiple
#' close-kinships
#'
#' @inheritParams find_ale_frqs
#' @inheritParams find_pss_frst_gt_prbs
#' @param cks Vector of codes representing kinships for which to calculate
#'   possible genopair probabilities
#'
#' @return Array with dimensions (number of possible genotypes x number of
#'   possible genotypes x number of loci x number of kinships)
#' @export
#'
#' @examples
find_pss_gp_prbs_kps = function(ind_gts, L, cks = c("HSP", "POP", "SP")) {
  # Find allele frequencies
  ale_frqs = find_ale_frqs(ind_gts)

  # Find possible genotype probabilities
  pss_gt_prbs = find_pss_gt_prbs(ale_frqs)

  # Find possible first genotype probabilities for genopairs
  pss_frst_gt_prbs = find_pss_frst_gt_prbs(pss_gt_prbs, L)

  # Find possible genopair probabilities for unrelated pairs
  pss_gp_prbs_ups = find_pss_gp_prbs_ups(pss_frst_gt_prbs)

  # Set genopair probabilities to NULL
  pss_gp_prbs_pops = pss_gp_prbs_sps = pss_gp_prbs_hsps = NULL

  # If parent-offspring pairs in kinship set selected
  if("POP" %in% cks) {
    # Find possible genopair probabilities for parent-offspring pairs
    pss_gp_prbs_pops = find_pss_gp_prbs_pops(ale_frqs, L, pss_frst_gt_prbs)
  }

  # If self-pairs in kinship set selected
  if("SP" %in% cks) {
    # Find possible genopair probabilities for self pairs
    pss_gp_prbs_sps = find_pss_gp_prbs_sps(pss_frst_gt_prbs, L)
  }

  # If half-sibling pairs in kinship set selected
  if("HSP" %in% cks) {
    # Find possible genopair probabilities for half-sibling pairs
    pss_gp_prbs_hsps = find_pss_gp_prbs_hsps(pss_gp_prbs_ups, pss_gp_prbs_pops)
  }

  # Set close kinships
  close_kps = c("HSP", "POP", "SP")

  # Combine genopair probabilities for selected kinships and return
  array(
    data = c(pss_gp_prbs_ups, pss_gp_prbs_hsps, pss_gp_prbs_pops,
             pss_gp_prbs_sps),
    dim = c(n_pss_gts, n_pss_gts, L, length(cks) + 1),
    dimnames = c(dimnames(pss_frst_gt_prbs),
                 list(kinship = c("UP", close_kps[close_kps %in% cks])))
  )
}

