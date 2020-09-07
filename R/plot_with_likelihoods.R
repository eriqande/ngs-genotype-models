

plot_with_likelihoods <- function(input, genos, rq, genolikes) {

  # and set, here, some values for their positions.
  bot <- 0.8  # distance below the ypos
  top <- 0.3  # distance below baseline
  w <- 0.23   # width of each bar
  full_h <- bot - top  # full height of a bar of likelihood 1
  genolikes2 <- genolikes %>%
    mutate(
      hg_f = factor(
        hyp_geno, levels = c("CC", "CT", "TT")
      ),
      w = w,
      full_h = full_h,
      xfl = xpos - 1.5 * w,  # x far left
      xfr = xpos + 1.5 * w,   # x far right
      ybot = ypos - bot,
      xl = xfl + (as.integer(hg_f) - 1.0) * w
    )

  count_sep = 0.06
  rd_summ <- rq %>%
    group_by(idx, sequenced_allele) %>%
    tally() %>%
    ungroup() %>%
    complete(
      idx = unique(genos$idx),
      sequenced_allele = c("C", "T"),
      fill = list(n = 0L)
    ) %>%
    left_join(
      genos %>%
        select(idx, xpos, ypos) %>%
        distinct(),
      by = "idx"
    ) %>%
    mutate(
      idx = as.integer(idx),
      sequenced_allele = as.character(sequenced_allele)
    ) %>%
    mutate(
      ypos = case_when(
        sequenced_allele == "C" ~ ypos - count_sep,
        TRUE ~ ypos - count_sep * 2
      )
    )

  ynum <- max(n_distinct(genos$ypos) - 1, 1)
  ylo <- min(genos$ypos) - 1.0
  yhi <- max(genos$ypos) + 0.2 / ynum

  geno_plot <- ggplot(genos, aes(x = xpos)) +
    geom_text(aes(label = geno, colour = geno, y = ypos)) +
    geom_text(aes(label = idx + 1, y = ypos + 0.08), size = 2.5) +
    geom_text(
      data = rd_summ,
      mapping = aes(label = n, y = ypos, colour = sequenced_allele)) +
    geom_rect(   # get the "bar plots" of the likelihoods
      data = genolikes2,
      aes(
        xmin = xl,
        xmax = xl + w,
        ymin = ybot,
        ymax = ybot + (scaled_likelihood * full_h),
        fill = hyp_geno
      )
    ) +
    geom_rect(  # get the dotted line perimeter on the bar plots
      data = genolikes2,
      aes(
        xmin = xfl,
        xmax = xfr,
        ymin = ybot,
        ymax = ybot + full_h
      ),
      linetype = "dotted",
      colour = "gray",
      fill = NA
    ) +
    geom_segment( # get the black line at the bottom
      data = genolikes2,
      aes(x = xfl, xend = xfr, y = ybot, yend = ybot),
      colour = "black"
    ) +
    ylim(ylo, yhi) +
    theme_void() +
    theme(
      legend.position = "none"
    ) +
    scale_colour_manual(
      values = c(
        CC = "blue",
        CT = "violet",
        TT = "red",
        C = "blue",
        T = "red"
      )
    ) +
    scale_fill_manual(
      values = c(
        CC = "blue",
        CT = "violet",
        TT = "red",
        C = "blue",
        T = "red"
      )
    ) +
    ggtitle("Sample of Genotypes:")

  numT <- sum(genos$g012)
  numC <- nrow(genos) * 2 - numT

  post_freq_plot <- ggplot() +
    stat_function(
      fun = dbeta,
      n = 1e4,
      args = list(
        shape1 = numT + input$alpha_1,
        shape2 = numC + input$alpha_2
      )
    ) +
    xlim(0,1) +
    xlab("p_T") +
    ylab("Posterior Density") +
    geom_vline(xintercept = input$p_T, colour = "red")

  # return the two plots so that we can modify them each as need be
  # later and combine as we want
  list(
    geno_plot = geno_plot,
    post_freq_plot = post_freq_plot,
    genolikes2 = genolikes2
  )

}
