# ============================================================================
# CALCULATE NBES
# ============================================================================

calculate_nbes <- function(
    high_df,
    low_df
) {
  
  high <- high_df %>%
    
    dplyr::select(
      
      mutation,
      
      log2E,
      log2E_norm
      
    ) %>%
    
    dplyr::rename(
      
      log2E_high =
        log2E,
      
      log2E_norm_high =
        log2E_norm
      
    )
  
  
  
  low <- low_df %>%
    
    dplyr::select(
      
      mutation,
      
      log2E,
      log2E_norm
      
    ) %>%
    
    dplyr::rename(
      
      log2E_low =
        log2E,
      
      log2E_norm_low =
        log2E_norm
      
    )
  
  
  
  nbes_df <- dplyr::inner_join(
    
    high,
    low,
    
    by = "mutation"
    
  )
  
  
  
  nbes_df <- nbes_df %>%
    
    dplyr::mutate(
      
      NBES =
        
        log2E_norm_high -
        
        log2E_norm_low
      
    )
  
  
  
  nbes_df <- nbes_df %>%
    
    dplyr::arrange(
      dplyr::desc(NBES)
    ) %>%
    
    dplyr::mutate(
      rank = dplyr::row_number()
    )
  
  
  
  summary <- list(
    
    n_high_variants =
      nrow(high),
    
    n_low_variants =
      nrow(low),
    
    n_variants_used =
      nrow(nbes_df)
    
  )
  
  
  
  list(
    
    nbes_df =
      nbes_df,
    
    summary =
      summary
    
  )
  
}

# ============================================================================
# CREATE NBES DATASET
# ============================================================================

create_nbes_dataset <- function(
    nbes_results,
    dataset_name,
    high_dataset,
    low_dataset
) {
  
  list(
    
    name = dataset_name,
    
    stage = "nbes",
    
    source = "nbes",
    
    created = Sys.time(),
    
    metadata = list(
      
      high_dataset =
        high_dataset,
      
      low_dataset =
        low_dataset,
      
      nbes_summary =
        nbes_results$summary
      
    ),
    
    results = list(
      
      nbes = nbes_results$nbes_df,
      
      position_summary =
        summarise_nbes_by_position(
          nbes_results$nbes_df
        )
      
    )
    
  )
  
}