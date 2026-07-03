# table for top-enriched variants
top_enriched_variants <- function(
    enrichment_df,
    n = 20
){
  
  enrichment_df %>%
    
    arrange(
      desc(log2E_norm)
    ) %>%
    
    slice_head(
      n = n
    )
  
}

# top depleted variants
top_depleted_variants <- function(
    enrichment_df,
    n = 20
){
  
  enrichment_df %>%
    
    arrange(
      log2E_norm
    ) %>%
    
    slice_head(
      n = n
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

# ============================================================================
# TOP POSITIVE NBES
# ============================================================================

top_positive_nbes <- function(
    nbes_df,
    n = 20
){
  
  nbes_df %>%
    
    arrange(
      desc(NBES)
    ) %>%
    
    slice_head(
      n = n
    )
  
}



# ============================================================================
# TOP NEGATIVE NBES
# ============================================================================

top_negative_nbes <- function(
    nbes_df,
    n = 20
){
  
  nbes_df %>%
    
    arrange(
      NBES
    ) %>%
    
    slice_head(
      n = n
    )
  
}

position_summary_enrichment <- function(dataset) {
  
  dataset$results$position_summary
  
}

position_summary_nbes <- function(dataset) {
  
  dataset$results$position_summary
  
}