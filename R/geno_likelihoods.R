
# Function to compute likelihoods for a single indiv
# and return them in a tibble.  Group on idx and geno to
# send stuff to this function
single_indiv_geno_like <- function(seq_err_prob, sequenced_allele) {

    scaled_likelihood <- numeric(3)
    scaled_likelihood[1] <- prod(
      (1 - seq_err_prob) * (sequenced_allele == "C") +
        (seq_err_prob * (sequenced_allele == "T"))
    )
    scaled_likelihood[2] <- 0.5 ^ length(sequenced_allele)
    scaled_likelihood[3] <- prod(
      (1 - seq_err_prob) * (sequenced_allele == "T") +
        (seq_err_prob * (sequenced_allele == "C"))
    )
    tibble(
      hyp_geno = c("CC", "CT", "TT"),  # this is the "hypothetical" genotype
      scaled_likelihood = scaled_likelihood / sum(scaled_likelihood)
    )
}


#' compute likelihoods for CC, CT, and TT from read depth data
#'
#' @param reads_and_quals a data frame of read depths
#' @param genotypes a data frame of the original genotypes.  This is here
#' in order to re-insert individuals that had now reads.
geno_likelihoods <- function(reads_and_quals, genotypes) {
  # get likelihoods for everyone with reads
  L1 <- reads_and_quals %>%
    group_by(idx, geno, ypos, xpos) %>%
    summarise(
      GLs = list(
        single_indiv_geno_like(seq_err_prob, sequenced_allele)
      )
    ) %>%
    unnest(GLs) %>%
    ungroup()

  # assign 1/3 likelihoods for everyone with no reads
  L2 <- genotypes %>%
    select(idx, geno, ypos, xpos) %>%
    anti_join(L1, by = "idx") %>%
    mutate(
      GLs = list(
        tibble(
          hyp_geno = c("CC", "CT", "TT"),
          scaled_likelihood = 1.0/3
        )
      )
    ) %>%
    unnest(GLs)

  bind_rows(
    L1,
    L2
  ) %>%
    arrange(idx, hyp_geno)
}
