compute_posteriors_mcmc <- function(
  input,
  genos,
  rq,
  genolikes
) {

  burn_in <- input$burn_sweeps[1]
  sweeps <- input$burn_sweeps[2]

  # first, get the genotype likelihoods in a matrix format:.
  # This is 3 rows and n columns
  GLmat <- genolikes %>%
    select(idx, hyp_geno, scaled_likelihood) %>%
    pivot_wider(
      names_from = hyp_geno,
      values_from = scaled_likelihood
    ) %>%
    select(-idx) %>%
    as.matrix() %>%
    t()

  # make a place to store p output
  p <- numeric(sweeps)


  # get a starting P-value by weighting things by the
  # likelihoods
  p[1] <- mean(GLmat * 0:2)

  # then get a starting set of posteriors from that
  gposts <- list()
  unnormo <- c(
    (1 - p[1]) ^ 2,
    2 * p[1] * (1 - p[1]),
    p[1] ^ 2
  ) * GLmat
  # normalize all those
  gposts[[1]] <- unnormo / rep(colSums(unnormo), each = 3)


  # now cycle over sweeps.  Unfortunately we can't use lapply, because
  # we are updating values, so we will just cycle over this in a for loop.

  for(t in 2:sweeps) {
    # first step is to simulate some zeds (i.e. latent genotypes)
    zeds <- apply(
      X = gposts[[t-1]],
      MARGIN = 2,
      FUN =  function(c) sample(x = 0:2, size = 1, replace = TRUE, prob = c)
    )
    # count up the number of 0 (C) alleles and the number of 1 (T) alleles
    num1 <- sum(zeds)
    num0 <- 2 * length(zeds) - num1

    # simulate a new p
    p[t] <- rbeta(1, input$alpha_1 + num1, input$alpha_2 + num0)

    # simulate a new set of posteriors
    unnormo <- c(
      (1 - p[t]) ^ 2,
      2 * p[t] * (1 - p[t]),
      p[t] ^ 2
    ) * GLmat
    # normalize all those
    gposts[[t]] <- unnormo / rep(colSums(unnormo), each = 3)
  }

  # now, get the mean posteriors for the genotypes
  X <- gposts[burn_in:sweeps]
  gp_array <- do.call(cbind, X)
  gp_array <- array(gp_array, dim=c(dim(X[[1]]), length(X)))
  GP <- colMeans(aperm(gp_array, c(3, 1, 2)), na.rm = TRUE)

  # turn those into a tibble, too.
  GP_tibble <- GP %>%
    t() %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(c("CC", "CT", "TT")) %>%
    mutate(idx = (1:n()) - 1L) %>%
    pivot_longer(
      names_to = "hyp_geno",
      values_to = "geno_posterior",
      -idx
    )

  # finally, make a tibble of the p-values
  p_tibble <- tibble(
    p = p
  ) %>%
    mutate(sweep = 1:n()) %>%
    select(sweep, p)

  # return both of those
  list(
    GP = GP_tibble,
    p = p_tibble
  )

}
