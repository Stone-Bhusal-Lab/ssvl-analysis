# ============================================================================
# COMPARISON MODULE - UI
# ============================================================================

mod_comparison_ui <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    # ------------------------------------------------------------------------
    # Dataset Selection
    # ------------------------------------------------------------------------
    
    fluidRow(
      
      column(
        6,
        
        selectInput(
          ns("dataset_a"),
          "Dataset A",
          choices = NULL
        )
      ),
      
      column(
        6,
        
        selectInput(
          ns("dataset_b"),
          "Dataset B",
          choices = NULL
        )
      )
      
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # QC Comparison
    # ------------------------------------------------------------------------
    
    h3("QC Comparison"),
    
    tableOutput(
      ns("qc_comparison")
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Variant Overlap
    # ------------------------------------------------------------------------
    
    h3("Variant Overlap"),
    
    tableOutput(
      ns("overlap_summary")
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Frequency Correlation
    # ------------------------------------------------------------------------
    
    h3("Comparison Statistics"),
    
    tableOutput(
      ns("comparison_stats")
    ),
    
    verbatimTextOutput(
      ns("jaccard")
    ),
    
    hr(),
    
    h3("Frequency Correlation"),
    
    plotOutput(
      ns("freq_plot"),
      height = "500px"
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Shared Variants
    # ------------------------------------------------------------------------
    
    h3("Shared Variants"),
    
    DT::DTOutput(
      ns("shared_variants")
    ),
    
    hr(),
    
    h3("Top Differential Variants"),
    
    DT::DTOutput(
      ns("differential_variants")
    )
    
  )
  
}
# ============================================================================
# COMPARISON MODULE - SERVER
# ============================================================================

mod_comparison_server <- function(
    id,
    datasets,
    active_dataset
){
  
  moduleServer(
    
    id,
    
    function(input, output, session){
      
      # ----------------------------------------------------------------------
      # Dataset Choices
      # ----------------------------------------------------------------------
      
      observe({
        
        choices <- names(datasets())
        
        updateSelectInput(
          session,
          "dataset_a",
          choices = choices
        )
        
        updateSelectInput(
          session,
          "dataset_b",
          choices = choices
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Dataset Reactives
      # ----------------------------------------------------------------------
      
      dataset_a <- reactive({
        
        req(input$dataset_a)
        
        datasets()[[input$dataset_a]]
        
      })
      
      dataset_b <- reactive({
        
        req(input$dataset_b)
        
        datasets()[[input$dataset_b]]
        
      })
      
      # ----------------------------------------------------------------------
      # QC Comparison
      # ----------------------------------------------------------------------
      
      output$qc_comparison <- renderTable({
        
        qa <- dataset_a()$results$qc
        qb <- dataset_b()$results$qc
        
        data.frame(
          
          Metric = c(
            "Total Reads",
            "Valid Inserts",
            "Extraction Rate (%)",
            "Filtered Haplotypes",
            "Single Mutants"
          ),
          
          Dataset_A = c(
            qa$total_reads,
            qa$valid_inserts,
            round(
              100 * qa$valid_inserts /
                qa$total_reads,
              1
            ),
            qa$unique_haplotypes_filtered,
            nrow(
              dataset_a()$results$single_mutants
            )
          ),
          
          Dataset_B = c(
            qb$total_reads,
            qb$valid_inserts,
            round(
              100 * qb$valid_inserts /
                qb$total_reads,
              1
            ),
            qb$unique_haplotypes_filtered,
            nrow(
              dataset_b()$results$single_mutants
            )
          )
          
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Shared Variant Table
      # ----------------------------------------------------------------------
      
      shared_variants <- reactive({
        
        a <- dataset_a()$results$single_mutants
        b <- dataset_b()$results$single_mutants
        
        inner_join(
          a,
          b,
          by = "mutation",
          suffix = c("_A", "_B")
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Overlap Summary
      # ----------------------------------------------------------------------
      
      output$overlap_summary <- renderTable({
        
        a <- dataset_a()$results$single_mutants
        b <- dataset_b()$results$single_mutants
        
        shared <- nrow(
          inner_join(
            a,
            b,
            by = "mutation"
          )
        )
        
        unique_a <- nrow(
          anti_join(
            a,
            b,
            by = "mutation"
          )
        )
        
        unique_b <- nrow(
          anti_join(
            b,
            a,
            by = "mutation"
          )
        )
        
        data.frame(
          
          Metric = c(
            "Shared Variants",
            "Unique To Dataset A",
            "Unique To Dataset B"
          ),
          
          Count = c(
            shared,
            unique_a,
            unique_b
          )
          
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Frequency Correlation
      # ----------------------------------------------------------------------
      
      freq_compare <- reactive({
        
        a <- dataset_a()$results$single_mutants %>%
          select(
            mutation,
            freq
          ) %>%
          rename(
            freq_a = freq
          )
        
        b <- dataset_b()$results$single_mutants %>%
          select(
            mutation,
            freq
          ) %>%
          rename(
            freq_b = freq
          )
        
        full_join(
          a,
          b,
          by = "mutation"
        ) %>%
          mutate(
            freq_a = replace_na(
              freq_a,
              0
            ),
            freq_b = replace_na(
              freq_b,
              0
            )
          )
        
      })
      
      # ----------------------------------------------------------------------
      # Comparison Statistics
      # ----------------------------------------------------------------------
      
      comparison_stats <- reactive({
        
        df <- freq_compare()
        
        overlapping <- df %>%
          filter(
            freq_a > 0,
            freq_b > 0
          )
        
        data.frame(
          
          Metric = c(
            "Shared Variants",
            "Pearson Correlation",
            "Spearman Correlation"
          ),
          
          Value = c(
            
            nrow(overlapping),
            
            round(
              cor(
                df$freq_a,
                df$freq_b,
                method = "pearson"
              ),
              3
            ),
            
            round(
              cor(
                df$freq_a,
                df$freq_b,
                method = "spearman"
              ),
              3
            )
            
          )
          
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Jaccard Similarity
      # ----------------------------------------------------------------------
      
      jaccard_similarity <- reactive({
        
        a <- dataset_a()$results$single_mutants$mutation
        
        b <- dataset_b()$results$single_mutants$mutation
        
        round(
          length(intersect(a, b)) /
            length(union(a, b)),
          3
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Differential Variants
      # ----------------------------------------------------------------------
      
      differential_variants <- reactive({
        
        df <- freq_compare()
        
        df %>%
          
          mutate(
            freq_difference =
              abs(freq_a - freq_b)
          ) %>%
          
          arrange(
            desc(freq_difference)
          )
        
      })
      
      # ----------------------------------------------------------------------
      # Comparison Statistics Table
      # ----------------------------------------------------------------------
      
      output$comparison_stats <- renderTable({
        
        comparison_stats()
        
      })
      # ----------------------------------------------------------------------
      # Jaccard Similarity
      # ----------------------------------------------------------------------
      
      output$jaccard <- renderText({
        
        paste(
          "Jaccard Similarity:",
          jaccard_similarity()
        )
        
      })
      
      output$freq_plot <- renderPlot({
        
        ggplot(
          
          freq_compare(),
          
          aes(
            freq_a,
            freq_b
          )
          
        ) +
          
          geom_point(
            alpha = 0.6
          ) +
          
          geom_abline(
            slope = 1,
            intercept = 0,
            colour = "red",
            linetype = "dashed"
          ) +
          
          theme_minimal() +
          
          labs(
            x = input$dataset_a,
            y = input$dataset_b,
            title = "Variant Frequency Correlation"
          )
        
      })
      
      # ----------------------------------------------------------------------
      # Shared Variant Table
      # ----------------------------------------------------------------------
      
      output$shared_variants <- DT::renderDT({
        
        shared_variants()
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # Differential Variants Table
      # ----------------------------------------------------------------------
      
      output$differential_variants <- DT::renderDT({
        
        differential_variants()
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
    }
    
  )
  
}