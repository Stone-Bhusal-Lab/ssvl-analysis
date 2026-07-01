# ============================================================================
# VISUALISATION MODULE - UI
# ============================================================================

mod_visualisation_ui <- function(id) {
  
  ns <- NS(id)
  
  fluidPage(
    
    # ------------------------------------------------------------------------
    # Dataset Information
    # ------------------------------------------------------------------------
    
    h3("Active Dataset"),
    
    verbatimTextOutput(
      ns("dataset_info")
    ),
    
    tableOutput(
      ns("qc")
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # QC Plots
    # ------------------------------------------------------------------------
    
    fluidRow(
      
      column(
        4,
        plotOutput(
          ns("haplotype_distribution")
        )
      ),
      
      column(
        4,
        plotOutput(
          ns("variant_classes")
        )
      ),
      
      column(
        4,
        plotOutput(
          ns("variant_reads")
        )
      )
      
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # QC Summary Table
    # ------------------------------------------------------------------------
    
    DT::DTOutput(
      ns("variant_qc")
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Variant Browser
    # ------------------------------------------------------------------------
    
    selectInput(
      ns("dataset_view"),
      "View",
      choices = c(
        "Single Mutants",
        "Filtered Variants",
        "Raw Variants",
        "Missense",
        "Indels",
        "Frameshifts"
      )
    ),
    
    DT::DTOutput(
      ns("variants")
    )
    
  )
  
}


# ============================================================================
# VISUALISATION MODULE - SERVER
# ============================================================================

mod_visualisation_server <- function(
    id,
    datasets,
    active_dataset
) {
  
  moduleServer(
    
    id,
    
    function(input, output, session) {
      
      # ----------------------------------------------------------------------
      # Active Dataset
      # ----------------------------------------------------------------------
      
      current_dataset <- reactive({
        
        req(active_dataset())
        
        datasets()[[active_dataset()]]
        
      })
      
      # ----------------------------------------------------------------------
      # Dataset Information
      # ----------------------------------------------------------------------
      
      output$dataset_info <- renderPrint({
        
        ds <- current_dataset()
        
        cat("Name:", ds$name, "\n")
        cat("Stage:", ds$stage, "\n")
        cat("Source:", ds$source, "\n")
        cat("Created:", ds$created, "\n")
        
        if (!is.null(ds$metadata)) {
          
          cat("\nAnalysis Parameters\n")
          cat("-------------------\n")
          
          if (!is.null(ds$metadata$min_count)) {
            cat("Min Count:", ds$metadata$min_count, "\n")
          }
          
          if (!is.null(ds$metadata$min_freq)) {
            cat("Min Frequency:", ds$metadata$min_freq, "\n")
          }
          
        }
        
      })
      
      # ----------------------------------------------------------------------
      # QC Metrics
      # ----------------------------------------------------------------------
      
      output$qc <- renderTable({
        
        ds <- current_dataset()
        
        as.data.frame(
          ds$results$qc
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Haplotype Distribution
      # ----------------------------------------------------------------------
      
      output$haplotype_distribution <- renderPlot({
        
        ds <- current_dataset()
        
        min_count <- DEFAULT_MIN_COUNT
        
        if (
          !is.null(ds$metadata) &&
          !is.null(ds$metadata$min_count)
        ) {
          min_count <- ds$metadata$min_count
        }
        
        min_count <- DEFAULT_MIN_COUNT
        
        if (
          !is.null(ds$metadata) &&
          !is.null(ds$metadata$min_count)
        ) {
          min_count <- ds$metadata$min_count
        }
        
        plot_haplotype_distribution(
          haplo_df_raw = ds$results$haplo_df_raw,
          min_count = min_count
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Variant Class Distribution
      # ----------------------------------------------------------------------
      
      output$variant_classes <- renderPlot({
        
        ds <- current_dataset()
        
        plot_variant_class_distribution(
          ds$results$haplo_df_raw
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Reads by Variant Class
      # ----------------------------------------------------------------------
      
      output$variant_reads <- renderPlot({
        
        ds <- current_dataset()
        
        plot_variant_class_reads(
          ds$results$haplo_df_raw
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Variant QC Summary
      # ----------------------------------------------------------------------
      
      output$variant_qc <- DT::renderDT({
        
        ds <- current_dataset()
        
        variant_qc_summary(
          ds$results$haplo_df_raw
        )
        
      },
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # Variant Table Selector
      # ----------------------------------------------------------------------
      
      variant_table <- reactive({
        
        ds <- current_dataset()
        
        switch(
          
          input$dataset_view,
          
          "Single Mutants" =
            ds$results$single_mutants,
          
          "Filtered Variants" =
            ds$results$haplo_df,
          
          "Raw Variants" =
            ds$results$haplo_df_raw,
          
          "Missense" =
            ds$results$missense,
          
          "Indels" =
            ds$results$indels,
          
          "Frameshifts" =
            ds$results$frameshifts,
          
          ds$results$single_mutants
          
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Variant Table
      # ----------------------------------------------------------------------
      
      output$variants <- DT::renderDT({
        
        variant_table()
        
      },
      options = list(
        pageLength = 25,
        scrollX = TRUE
      ))
      
    }
    
  )
  
}