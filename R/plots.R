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
      title = "Haplotype Read Count Distribution (pre-filter)",
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

# enirchment distribution
plot_enrichment_distribution <- function(
    enrichment_df
){
  
  ggplot(
    enrichment_df,
    aes(log2E_norm)
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    
    theme_minimal() +
    
    labs(
      title = "WT-Normalised Enrichment",
      x = "log2E_norm",
      y = "Number of Variants"
    )
  
}
# ============================================================================
# POSITION ENRICHMENT PROFILE
# ============================================================================

plot_position_enrichment <- function(
    position_summary
){
  
  ggplot(
    position_summary,
    aes(
      x = position,
      y = mean_log2E_norm
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_line() +
    
    geom_point() +
    
    theme_minimal() +
    
    labs(
      title = "Position Enrichment Profile",
      x = "Position",
      y = "Mean log2E_norm"
    )
  
}
# ============================================================================
# NBES DISTRIBUTION
# ============================================================================

plot_nbes_distribution <- function(
    nbes_df
){
  
  ggplot(
    nbes_df,
    aes(NBES)
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    theme_minimal() +
    
    labs(
      title = "NBES Distribution",
      x = "NBES",
      y = "Count"
    )
  
}



# ============================================================================
# POSITION NBES PROFILE
# ============================================================================

plot_position_nbes <- function(
    position_summary
){
  
  ggplot(
    position_summary,
    aes(
      x = position,
      y = mean_nbes
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_line() +
    
    geom_point() +
    
    theme_minimal() +
    
    labs(
      title = "Position NBES Profile",
      x = "Position",
      y = "Mean NBES"
    )
  
}

# ============================================================================
# enrichment heatmap
# ============================================================================

make_enrichment_heatmap <- function(
    df,
    value_col,
    title,
    ref_protein,
    reverse = FALSE
) {
  
  # ===============================
  # STEP 1: Extract + collapse ✅
  # ===============================
  heat <- df %>%
    filter(mutation != "WT") %>%
    mutate(
      position = as.numeric(gsub("[^0-9]", "", mutation)),
      alt = gsub(".*[0-9]+", "", mutation),
      value = .data[[value_col]]
    ) %>%
    group_by(position, alt) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # ===============================
  # STEP 2: Full AA space ✅
  # ===============================
  aa_levels <- c(
    "A","R","N","D","C","Q","E","G","H","I",
    "L","K","M","F","P","S","T","W","Y","V","del"
  )
  
  full_grid <- expand.grid(
    position = 1:nchar(DEFAULT_REF_PROTEIN),
    alt = aa_levels
  )
  
  # ===============================
  # STEP 3: Join ✅
  # ===============================
  heat <- full_grid %>%
    left_join(heat, by = c("position", "alt"))
  
  # leave missing as NA
  heat$value[is.na(heat$value)] <- NA
  
  # ===============================
  # STEP 4: Wide format ✅
  # ===============================
  heat_wide <- heat %>%
    pivot_wider(
      names_from = alt,
      values_from = value
    )
  
  # ===============================
  # STEP 5: Build matrix ✅
  # ===============================
  mat <- as.data.frame(heat_wide)
  rownames(mat) <- mat$position
  mat <- mat[, -1]
  
  mat_t <- t(as.matrix(mat))   # ✅ critical transpose
  
  # ===============================
  # STEP 6: Order columns ✅
  # ===============================
  mat_plot <- mat_t[, order(as.numeric(colnames(mat_t)))]
  
  # ===============================
  # STEP 7: Order rows ✅
  # ===============================
  rows <- rownames(mat_plot)
  aa_rows <- sort(setdiff(rows, "del"))
  row_order <- c(aa_rows, "del"[ "del" %in% rows ])
  mat_plot <- mat_plot[row_order, ]
  
  # ===============================
  # STEP 8: WT labels ✅
  # ===============================
  wt_aa <- strsplit(DEFAULT_REF_PROTEIN, "")[[1]]
  wt_labels <- wt_aa[as.numeric(colnames(mat_plot))]
  
  top_anno <- HeatmapAnnotation(
    WT = anno_text(
      wt_labels,
      rot = 0,
      gp = gpar(fontsize = 10),
      just = "center"
    ),
    height = unit(0.1, "cm")
  )
  
  # ===============================
  # STEP 9: Colour scale ✅
  # ===============================
  max_abs <- max(abs(mat_plot), na.rm = TRUE)
  if (!is.finite(max_abs)) max_abs <- 1
  
  if (!reverse) {
    col_fun <- colorRamp2(
      c(-max_abs, 0, max_abs),
      c("#2166AC", "#FFFFFF", "#B2182B")   # blue → white → red
    )
  } else {
    col_fun <- colorRamp2(
      c(-max_abs, 0, max_abs),
      c("#B2182B", "#FFFFFF", "#2166AC")   # ✅ reversed
    )
  }
  Heatmap(
    mat_plot,
    name = value_col,
    col = col_fun,
    na_col = "grey40",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_centered = TRUE,
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
    top_annotation = top_anno,
    column_title = paste0(title, "\n"),
    column_title_gp = gpar(fontsize = 15)
  )
  
}

# ============================================================================
# PLOT ENRICHMENT HEATMAP
# ============================================================================

plot_enrichment_heatmap <- function(
    enrichment_df
){
  
  make_enrichment_heatmap(
    df = enrichment_df,
    value_col = "log2E_norm",
    title = "WT-Normalised Enrichment Heatmap",
    reverse = FALSE
  )
  
}


# ============================================================================
# NBES HEATMAP
# ============================================================================

plot_nbes_heatmap <- function(
    nbes_df,
    ref_protein
){
  
  make_enrichment_heatmap(
    df = nbes_df,
    value_col = "NBES",
    title = "NBES Heatmap",
    ref_protein = ref_protein
  )
  
}