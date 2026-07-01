plot_haplotype_distribution <- function(haplo_df_raw, min_count) {
  
  count_dist <- haplo_df_raw %>%
    dplyr::count(count)
  
  ggplot2::ggplot(count_dist, ggplot2::aes(x = count, y = n)) +
    ggplot2::geom_point(
      alpha = 0.4,
      size = 1.5,
      colour = "steelblue"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(
      xintercept = min_count,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Read Count",
      y = "Number of Haplotypes",
      title = "Read Count Distribution"
    )
}