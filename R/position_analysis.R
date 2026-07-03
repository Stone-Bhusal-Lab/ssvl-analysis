# ============================================================
# position-level summary functions
# ============================================================
# ------------------------------------------------------------
# Extract position position from mutation string
# Examples:
# WT     -> NA
# A23V   -> 23
# R102K  -> 102
# G56*   -> 56
# ------------------------------------------------------------

extract_position_position <- function(mutation) {
  
  position <- gsub("[A-Za-z\\*]", "", mutation)
  
  suppressWarnings(as.numeric(position))
}


# ------------------------------------------------------------
# Add position column to variant table
# ------------------------------------------------------------

add_position_column <- function(df,
                               mutation_col = "mutation") {
  
  if (!mutation_col %in% names(df)) {
    stop("Column '", mutation_col, "' not found.")
  }
  
  df %>%
    mutate(
      position = ifelse(
        .data[[mutation_col]] == "WT",
        NA_real_,
        extract_position_position(.data[[mutation_col]])
      )
    )
}


# ------------------------------------------------------------
# Enrichment position summary
#
# Requires:
# mutation
# log2E
#
# Optional:
# log2E_norm
# ------------------------------------------------------------

summarise_enrichment_by_position <- function(df) {
  
  required_cols <- c("mutation", "log2E")
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df <- add_position_column(df)
  
  df %>%
    filter(mutation != "WT") %>%
    group_by(position) %>%
    summarise(
      n_variants = n(),
      
      mean_log2E = mean(log2E, na.rm = TRUE),
      median_log2E = median(log2E, na.rm = TRUE),
      sd_log2E = sd(log2E, na.rm = TRUE),
      
      min_log2E = min(log2E, na.rm = TRUE),
      max_log2E = max(log2E, na.rm = TRUE),
      
      mean_abs_log2E = mean(abs(log2E), na.rm = TRUE),
      max_abs_log2E = max(abs(log2E), na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    arrange(position)
}


# ------------------------------------------------------------
# Enrichment summary using WT-normalised values
#
# Requires:
# mutation
# log2E_norm
# ------------------------------------------------------------

summarise_enrichment_norm_by_position <- function(df) {
  
  required_cols <- c("mutation", "log2E_norm")
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df <- add_position_column(df)
  
  df %>%
    filter(mutation != "WT") %>%
    group_by(position) %>%
    summarise(
      n_variants = n(),
      
      mean_log2E_norm = mean(log2E_norm, na.rm = TRUE),
      median_log2E_norm = median(log2E_norm, na.rm = TRUE),
      sd_log2E_norm = sd(log2E_norm, na.rm = TRUE),
      
      min_log2E_norm = min(log2E_norm, na.rm = TRUE),
      max_log2E_norm = max(log2E_norm, na.rm = TRUE),
      
      mean_abs_log2E_norm = mean(abs(log2E_norm), na.rm = TRUE),
      max_abs_log2E_norm = max(abs(log2E_norm), na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    arrange(position)
}


# ------------------------------------------------------------
# NBES position summary
#
# Requires:
# mutation
# nbes
# ------------------------------------------------------------

summarise_nbes_by_position <- function(df) {
  
  required_cols <- c("mutation", "NBES")
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df <- add_position_column(df)
  
  df %>%
    filter(mutation != "WT") %>%
    group_by(position) %>%
    summarise(
      n_variants = n(),
      
      mean_nbes = mean(NBES, na.rm = TRUE),
      median_nbes = median(NBES, na.rm = TRUE),
      sd_nbes = sd(NBES, na.rm = TRUE),
      
      min_nbes = min(NBES, na.rm = TRUE),
      max_nbes = max(NBES, na.rm = TRUE),
      
      mean_abs_nbes = mean(abs(NBES), na.rm = TRUE),
      max_abs_nbes = max(abs(NBES), na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    arrange(position)
}