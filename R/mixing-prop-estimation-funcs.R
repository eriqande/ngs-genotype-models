

#### Start off defining all of our helper functions ####

# for plotting the allele frequencies
allele_freq_scatters <-  function(p1, p2) {
  tibble(p1 = p1, p2 = p2) %>%
    ggplot(aes(x = p1, y = p2)) +
    geom_point(col = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted")
}


# compute likelihoods for individuals from different types of genotype matrix input,
# as we will use the genotype likelihoods or maximum likelihood genotypes here as well
ret_lik <- function(x, which = 1, GP) {
  logl1 <- sum(log(colSums(x * GP[[1]])))
  logl2 <- sum(log(colSums(x * GP[[2]])))
  lmax <- pmax(logl1, logl2)
  lik1 <- exp(logl1 - lmax) / (exp(logl1 - lmax) + exp(logl2 - lmax))
  lik2 <- exp(logl2 - lmax) / (exp(logl1 - lmax) + exp(logl2 - lmax))
  if(which == 1) {
    return(lik1)
  }
  if(which == 2) {
    return(lik2)
  }
  else {stop("which must be one or two")}
}



# define a function for getting the posterior for the mixing
# proportion from the L1 and L2, which are the likelihoods that you want to use.
mix_post_func <- function(L1, L2, pts = seq(0.01, 0.99, by = 0.005)) {
  mixes <- pts
  mix_post <- unlist(lapply(mixes, function(m) {
    sum(log(m * L1 + (1 - m) * L2))
  }))
  mix_post <- exp(mix_post - max(mix_post))

  # now we need to normalize this so it integrates to 1 over these little segments
  (length(pts) - 1) * mix_post / sum(mix_post)
}


#' Given read depths in a two-vector (one value for each allele),
#' this function returns the genotype likelihoods
#' @param R a two vector: number of reads of allele 1 and number of allele 2.
#' @return A three vector: the likelihood that the true genotype is 11, 12, 22 where 1
#' is the allele for which the allele frequency is given.
glike_mix <- function(R) {
  n <- sum(R)
  l11 <- dbinom(R[1], n, prob = 1)
  l12 <- dbinom(R[1], n, prob = 0.5)
  l22 <- dbinom(R[1], n, prob = 0)

  c(l11, l12, l22) / (l11 + l12 + l22)
}


# here is a function to look at a genotype columns (a three vector with one 1 in it),
# and then simulate reads from that in a two vector, and then calculate the
# genotype likeihood of that with glike_mix. R is the average read depth
sim_reads_etc <- function(x, R) {
  RD = rpois(1, R)
  if(x[1] == 1) {
    gr <- c(RD, 0)
  } else if(x[3] == 1) {
    gr <- c(0, RD)
  } else {
    y <- rbinom(1, RD, 0.5)
    gr <- c(y, RD - y)
  }
  glike_mix(gr)
}



# here is another function to return a three-vector with a 1 at the maximum
# likelihood position (or a 0 if it is all missing)
max_geno_from_likes <- function(x) {
  ret <- c(0, 0, 0)
  if(x[1] == x[2] & x[2] == x[3]) {
    return(x)  # in this case, we say that any genotype is equally likely.  This will work out this way
  }
  z <- which.max(x)
  ret[z] <- 1
  ret
}


#### Now, here is a function to return a list of plots and things ####

simulate_data_and_inferences <- function(
  N = 500,
  L = 100,
  f = 0.1,
  mix = 0.9,
  AveRD = 0.1
) {

  # simulate allele frequencies
  p1 <- rbeta(L, 1.5, 1.5)
  p2 <- rbeta(L, p1 * (1 - f) / f, (1 - p1) * (1 - f) / f)

  # make a plot of the allele frequencies
  afreq_plot <- tibble(p1 = p1, p2 = p2) %>%
    ggplot(aes(x = p1, y = p2)) +
    geom_point(col = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    ggtitle("Allele Frequencies (Pop 1 on x; Pop 2 on y)")

  # get the genotype freqs in some matrices
  GP1 <- matrix(
    c(p1 ^ 2, 2 * p1 * (1 - p1), (1 - p1) ^ 2),
    byrow = TRUE,
    nrow = 3
  )
  GP2 <- matrix(
    c(p2 ^ 2, 2 * p2 * (1 - p2), (1 - p2) ^ 2),
    byrow = TRUE,
    nrow = 3
  )
  GP <- list(GP1, GP2)

  # Simulate the population memberships
  Z <- sample(1:2, size = N, replace = TRUE, prob = c(mix, 1 - mix))

  # Simulated the genotypes as matrices
  GM <- lapply(1:N, function(i) {
    gp <- GP[[Z[i]]]
    apply(gp, 2, function(x) rmultinom(1, 1, x))
  })

  # calculate the genotype likelihoods
  samples <- tibble(
    idx = 1:N,
    Z = Z,
    GM = GM
  ) %>%
    mutate(
      lik1 = map_dbl(.x = GM, ret_lik, which = 1, GP = GP),
      lik2 = map_dbl(.x = GM, ret_lik, which = 2, GP = GP)
    )

  # Simulate read depths and compute the genotype likelihoods from those
  # and also identify the maximum likelihood genotype
  samp2 <- samples %>%
    mutate(
      read_likes = map(
        .x = GM,
        .f = function(x) apply(x, 2, sim_reads_etc, R = AveRD)
      ),
      # while we are at it, let's just assign the maximum likelihood genotype,
      # and also the genotype posteriors given that they are from each of the
      # the two populations in turn
      max_like_geno = map(
        .x = read_likes,
        .f = function(x) apply(x, 2, max_geno_from_likes)
      )
    )


  # calculate genotypes by three different methods
  samp3 <- samp2 %>%
    mutate(
      like1_from_max_like_geno = map_dbl(.x = max_like_geno, ret_lik, which = 1, GP = GP),
      like2_from_max_like_geno = map_dbl(.x = max_like_geno, ret_lik, which = 2, GP = GP),
      like1_from_read_likes = map_dbl(.x = read_likes, ret_lik, which = 1, GP = GP),
      like2_from_read_likes = map_dbl(.x = read_likes, ret_lik, which = 2, GP = GP)
    )

  # compute the posterior for the mixing proportion
  pts <- seq(0, 1, by = 0.005)
  mix_posts <- tibble(
    mix_prop = pts,
    from_actual_geno = mix_post_func(samp3$lik1, samp3$lik2, pts),
    from_ml_geno = mix_post_func(samp3$like1_from_max_like_geno, samp3$like2_from_max_like_geno, pts),
    from_geno_likes = mix_post_func(samp3$like1_from_read_likes, samp3$like2_from_read_likes, pts)
  )

  # then make a plot of the mixing proportion posteriors
  mixing_proportion_posteriors <- mix_posts %>%
    rename(
      from_the_actual_genotypes = from_actual_geno,
      from_max_likelihood_assigned_genotypes = from_ml_geno,
      from_genotype_likelihoods = from_geno_likes
    ) %>%
    pivot_longer(cols = starts_with("from"), values_to = "posterior", names_to = "method") %>%
    ggplot(aes(x = mix_prop, y = posterior, colour = method)) +
    geom_line(size = 1.0) +
    geom_vline(xintercept = mix, linetype = "dashed") +
    ggtitle("Mixing Proportion Posteriors from Three Methods (True Value at Dashed Vertical Line)")


  # Finally, calculate the posterior probs of origin for each fish
  pop_1_assign_probs <- samp3 %>%
    mutate(
      from_the_actual_genotypes = map2_dbl(
        .x = lik1,
        .y = lik2,
        .f = function(x, y) {
          sum(mix_posts$from_actual_geno * (mix_posts$mix_prop * x) / (mix_posts$mix_prop * x + (1 - mix_posts$mix_prop) * y)) /
            sum(mix_posts$from_actual_geno)
        }
      ),
      from_max_likelihood_assigned_genotypes = map2_dbl(
        .x = like1_from_max_like_geno,
        .y = like2_from_max_like_geno,
        .f = function(x, y) {
          sum(mix_posts$from_ml_geno * (mix_posts$mix_prop * x) / (mix_posts$mix_prop * x + (1 - mix_posts$mix_prop) * y)) /
            sum(mix_posts$from_ml_geno)
        }
      ),
      from_genotype_likelihoods = map2_dbl(
        .x = like1_from_read_likes,
        .y = like2_from_read_likes,
        .f = function(x, y) {
          sum(mix_posts$from_geno_likes * (mix_posts$mix_prop * x) / (mix_posts$mix_prop * x + (1 - mix_posts$mix_prop) * y)) /
            sum(mix_posts$from_geno_likes)
        }
      )
    ) %>%
    select(idx, Z, starts_with("from")) %>%
    rename(true_population_of_origin = Z) %>%
    mutate(true_population_of_origin = str_c("Pop_", true_population_of_origin))

  # now, make a plot of those assignment probs
  ass_prob_histos <- pop_1_assign_probs %>%
    pivot_longer(
      cols = starts_with("from"),
      names_to = "method",
      values_to = "Pop1_posterior_prob"
      ) %>%
    ggplot(aes(x = Pop1_posterior_prob, fill = true_population_of_origin)) +
    geom_histogram(bins = 30) +
    facet_grid(true_population_of_origin ~ method, scales = "free_y") +
    ggtitle("Posterior probabilities of membership from Population 1")


  # finally, return some plots
  list(
    afreq_plot = afreq_plot,
    mixing_proportion_posteriors = mixing_proportion_posteriors,
    ass_prob_histos = ass_prob_histos
  )

}
