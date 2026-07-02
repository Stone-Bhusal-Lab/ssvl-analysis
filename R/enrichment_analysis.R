# ============================================================================
# CALCULATE ENRICHMENT
# ============================================================================

calculate_enrichment <- function(
    presort_df,
    selected_df
) {
  
  presort <- presort_df %>%
    dplyr::select(
      mutation,
      freq
    ) %>%
    dplyr::rename(
      freq_presort = freq
    )
  
  selected <- selected_df %>%
    dplyr::select(
      mutation,
      freq
    ) %>%
    dplyr::rename(
      freq_selected = freq
    )
  
  # --------------------------------------------------------------------------
  # Variants absent from pre-sort
  # --------------------------------------------------------------------------
  
  dropped_variants <- dplyr::anti_join(
    selected,
    presort,
    by = "mutation"
  )
  
  # --------------------------------------------------------------------------
  # Only keep variants present in pre-sort library
  # --------------------------------------------------------------------------
  
  enrichment_df <- dplyr::inner_join(
    presort,
    selected,
    by = "mutation"
  )
  
  # --------------------------------------------------------------------------
  # Dynamic pseudocount
  # --------------------------------------------------------------------------
  
  all_freqs <- c(
    enrichment_df$freq_presort,
    enrichment_df$freq_selected
  )
  
  all_freqs <- all_freqs[
    all_freqs > 0
  ]
  
  pseudocount <- min(all_freqs) / 10
  
  # --------------------------------------------------------------------------
  # Calculate raw enrichment
  # --------------------------------------------------------------------------
  
  enrichment_df <- enrichment_df %>%
    
    dplyr::mutate(
      
      log2E = log2(
        
        (freq_selected + pseudocount) /
          
          (freq_presort + pseudocount)
        
      )
      
    )
  
  # --------------------------------------------------------------------------
  # WT normalization
  # --------------------------------------------------------------------------
  
  wt_log2E <- enrichment_df %>%
    
    dplyr::filter(
      mutation == "WT"
    ) %>%
    
    dplyr::pull(
      log2E
    )
  
  if(length(wt_log2E) != 1){
    
    stop(
      "WT variant not found in enrichment dataset."
    )
    
  }
  
  enrichment_df <- enrichment_df %>%
    
    dplyr::mutate(
      
      log2E_norm =
        
        log2E - wt_log2E
      
    )
  
  # --------------------------------------------------------------------------
  # Ranking
  # --------------------------------------------------------------------------
  
  enrichment_df <- enrichment_df %>%
    
    dplyr::arrange(
      dplyr::desc(
        log2E_norm
      )
    ) %>%
    
    dplyr::mutate(
      rank = dplyr::row_number()
    )
  
  # --------------------------------------------------------------------------
  # Summary
  # --------------------------------------------------------------------------
  
  summary <- list(
    
    n_presort_variants =
      nrow(presort),
    
    n_selected_variants =
      nrow(selected),
    
    n_variants_used =
      nrow(enrichment_df),
    
    n_variants_dropped =
      nrow(dropped_variants),
    
    pct_variants_dropped =
      round(
        100 *
          nrow(dropped_variants) /
          max(1, nrow(selected)),
        2
      ),
    
    wt_log2E =
      wt_log2E,
    
    pseudocount =
      pseudocount
    
  )
  
  list(
    
    enrichment_df =
      enrichment_df,
    
    summary =
      summary,
    
    dropped_variants =
      dropped_variants
    
  )
  
}



# ============================================================================
# CREATE ENRICHMENT DATASET
# ============================================================================

create_enrichment_dataset <- function(
    enrichment_results,
    dataset_name,
    presort_name,
    selected_name
) {
  
  list(
    
    name = dataset_name,
    
    stage = "enrichment",
    
    source = "enrichment",
    
    created = Sys.time(),
    
    metadata = list(
      
      presort_dataset =
        presort_name,
      
      selected_dataset =
        selected_name,
      
      wt_log2E =
        enrichment_results$summary$wt_log2E,
      
      pseudocount =
        enrichment_results$summary$pseudocount,
      
      enrichment_summary =
        enrichment_results$summary
      
    ),
    
    results = list(
      
      enrichment =
        enrichment_results$enrichment_df,
      
      dropped_variants =
        enrichment_results$dropped_variants
      
    )
    
  )
  
}