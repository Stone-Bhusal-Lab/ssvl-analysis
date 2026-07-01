# Haplotype distribution pre-filter
plot_haplotype_distribution <- function(
    haplo_df_raw,
    min_count = NULL
) {
  
  if (is.null(haplo_df_raw) || nrow(haplo_df_raw) == 0) {
    
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle("No haplotypes available")
    )
    
  }
  
  if (!"count" %in% names(haplo_df_raw)) {
    
    stop("haplo_df_raw must contain a 'count' column")
    
  }
  
  count_dist <- haplo_df_raw %>%
    dplyr::filter(count > 0) %>%
    dplyr::count(count, name = "n_haplotypes")
  
  p <- ggplot2::ggplot(
    count_dist,
    ggplot2::aes(
      x = count,
      y = n_haplotypes
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.5,
      size = 1.5,
      colour = "steelblue"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Haplotype Read Count Distribution",
      x = "Read Count",
      y = "Number of Haplotypes"
    )
  
  if (!is.null(min_count)) {
    
    p <- p +
      ggplot2::geom_vline(
        xintercept = min_count,
        colour = "red",
        linetype = "dashed",
        linewidth = 0.8
      )
    
  }
  
  p
}
# Plot distribution of variant classes
plot_variant_class_distribution <- function(haplo_df) {
  
  if (is.null(haplo_df) || nrow(haplo_df) == 0) {
    
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle("No variant data available")
    )
    
  }
  
  if (!"variant_class" %in% names(haplo_df)) {
    
    stop("haplo_df must contain a 'variant_class' column")
    
  }
  
  class_counts <- haplo_df %>%
    dplyr::count(
      variant_class,
      name = "n_variants"
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_variants)
    )
  
  ggplot2::ggplot(
    class_counts,
    ggplot2::aes(
      x = reorder(variant_class, n_variants),
      y = n_variants,
      fill = variant_class
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Variant Class Distribution",
      x = NULL,
      y = "Number of Variants"
    )
}
# Read-Supported Composition
plot_variant_class_reads <- function(haplo_df) {
  
  class_counts <- haplo_df %>%
    dplyr::group_by(variant_class) %>%
    dplyr::summarise(
      total_reads = sum(count),
      .groups = "drop"
    )
  
  ggplot2::ggplot(
    class_counts,
    ggplot2::aes(
      x = reorder(variant_class, total_reads),
      y = total_reads,
      fill = variant_class
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Reads by Variant Class",
      x = NULL,
      y = "Total Reads"
    )
}
# Variant QC Summary
variant_qc_summary <- function(haplo_df) {
  
  total_reads <- sum(haplo_df$count)
  
  haplo_df %>%
    dplyr::group_by(variant_class) %>%
    dplyr::summarise(
      n_variants = dplyr::n(),
      total_reads = sum(count),
      pct_reads = 100 * sum(count) / total_reads,
      mean_frequency = mean(freq),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(total_reads)
    )
  
}