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
    
    hr(),
    
    h4("Plot Export"),
    
    selectInput(
      ns("plot_to_export"),
      "Plot",
      choices = c(
        "Position NBES",
        "Position Enrichment"
      )
    ),
    
    textInput(
      ns("export_title"),
      "Plot Title",
      value = ""
    ),
    
    textInput(
      ns("export_x"),
      "X Axis Title",
      value = ""
    ),
    
    textInput(
      ns("export_y"),
      "Y Axis Title",
      value = ""
    ),
    
    numericInput(
      ns("export_width"),
      "Width (cm)",
      value = 20,
      min = 1
    ),
    
    numericInput(
      ns("export_height"),
      "Height (cm)",
      value = 15,
      min = 1
    ),
    
    
    numericInput(
      ns("export_dpi"),
      "DPI",
      value = 300,
      min = 72
    ),
    
    selectInput(
      ns("export_format"),
      "Format",
      choices = c(
        "png",
        "pdf",
        "svg"
      )
    ),
    
    downloadButton(
      ns("export_plot"),
      "Export Plot"
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Mutation Dataset QC
    # ------------------------------------------------------------------------
    
    conditionalPanel(
      
      condition = sprintf(
        "output['%s'] == 'mutation'",
        ns("dataset_stage")
      ),
      
      tableOutput(
        ns("qc")
      ),
      
      hr(),
      
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
      
      DT::DTOutput(
        ns("variant_qc")
      ),
      
      hr()
      
    ),
    
    # ------------------------------------------------------------------------
    # Enrichment Visualisation
    # ------------------------------------------------------------------------
    
    conditionalPanel(
      
      condition = sprintf(
        "output['%s'] == 'enrichment'",
        ns("dataset_stage")
      ),
      
      h3("Enrichment Summary"),
      
      verbatimTextOutput(
        ns("enrichment_summary")
      ),
      
      fluidRow(
        
        column(
          6,
          plotOutput(
            ns("enrichment_distribution")
          )
        ),
        
        column(
          6,
          plotOutput(
            ns("position_enrichment")
          )
        )
        
      ),
      
      hr(),
      
      h4("Enrichment Heatmap"),
      
      plotOutput(
        ns("enrichment_heatmap"),
        height = "700px"
      ),
      
      hr(),
      
      hr(),
      
      h4("Top Enriched Variants"),
      
      DT::DTOutput(
        ns("top_enriched")
      ),
      
      hr(),
      
      h4("Top Depleted Variants"),
      
      DT::DTOutput(
        ns("top_depleted")
      ),
      
      hr(),
      
      h4("Position Summary"),
      
      DT::DTOutput(
        ns("position_enrichment_summary")
      )
      
    ),
    
    # ------------------------------------------------------------------------
    # NBES Visualisation
    # ------------------------------------------------------------------------
    
    conditionalPanel(
      
      condition = sprintf(
        "output['%s'] == 'nbes'",
        ns("dataset_stage")
      ),
      
      h3("NBES Summary"),
      
      verbatimTextOutput(
        ns("nbes_summary")
      ),
      
      fluidRow(
        
        column(
          6,
          plotOutput(
            ns("nbes_distribution")
          )
        ),
        
        column(
          6,
          plotOutput(
            ns("position_nbes")
          )
        )
        
      ),
      
      hr(),
      
      h4("NBES Heatmap"),
      
      plotOutput(
        ns("nbes_heatmap"),
        height = "700px"
      ),
      
      hr(),
      
      h4("Position Drilldown"),
      
      numericInput(
        ns("selected_position"),
        "Position",
        value = 1,
        min = 1
      ),
      
      DT::DTOutput(
        ns("position_variants")
      ),
      
      hr(),
      
      h4("Top Positive NBES"),
      
      DT::DTOutput(
        ns("top_positive_nbes")
      ),
      
      hr(),
      
      h4("Top Negative NBES"),
      
      DT::DTOutput(
        ns("top_negative_nbes")
      ),
      
      hr(),
      
      h4("Position Summary"),
      
      DT::DTOutput(
        ns("position_nbes_summary")
      ),
      
      
      
    ),
    
    # ------------------------------------------------------------------------
    # Variant Browser
    # ------------------------------------------------------------------------
    
    conditionalPanel(
      
      condition = sprintf(
        "output['%s'] == 'mutation'",
        ns("dataset_stage")
      ),
      
      selectInput(
        ns("dataset_view"),
        "View",
        choices = c(
          "Single Mutants",
          "Filtered Variants",
          "Missense",
          "Indels",
          "Frameshifts"
        )
      )
      
    ),
    
    conditionalPanel(
      
      condition = sprintf(
        "output['%s'] != 'mutation'",
        ns("dataset_stage")
      ),
      
      hr(),
      
      h4("Full Dataset"),
      
      DT::DTOutput(
        ns("variants")
      )
      
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
      
      dataset_stage <- reactive({
        
        current_dataset()$stage
        
      })
      
      enrichment_df <- reactive({
        
        ds <- current_dataset()
        
        req(ds$stage == "enrichment")
        
        ds$results$enrichment
        
      })
      
      nbes_df <- reactive({
        
        ds <- current_dataset()
        
        req(ds$stage == "nbes")
        
        ds$results$nbes
        
      })
      
      position_summary <- reactive({
        
        ds <- current_dataset()
        
        req(
          ds$stage %in% c(
            "enrichment",
            "nbes"
          )
        )
        
        ds$results$position_summary
        
      })
      
      export_plot_object <- reactive({
        
        ds <- current_dataset()
        
        req(
          ds$stage %in% c(
            "enrichment",
            "nbes"
          )
        )
        
        if (
          input$plot_to_export == "Position NBES" &&
          ds$stage == "nbes"
        ) {
          
          return(
            plot_position_nbes(
              position_summary()
            )
          )
          
        }
        
        if (
          input$plot_to_export == "Position Enrichment" &&
          ds$stage == "enrichment"
        ) {
          
          return(
            plot_position_enrichment(
              position_summary()
            )
          )
          
        }
        
        NULL
        
      })
      
      output$export_plot <- downloadHandler(
        
        filename = function() {
          paste0(
            "plot_",
            Sys.Date(),
            ".",
            input$export_format
          )
        },
        
        content = function(file) {
          
          p <- export_plot_object()
          
          req(p)
          
          p <- p +
            
            labs(
              title = input$export_title,
              x = input$export_x,
              y = input$export_y
            )
          
          ggsave(
            filename = file,
            plot = p,
            width = input$export_width,
            height = input$export_height,
            units = "cm",
            dpi = input$export_dpi
          )
          
        }
        
      )
      
      position_variants <- reactive({
        
        req(
          nbes_df(),
          input$selected_position
        )
        
        nbes_df() %>%
          
          mutate(
            position = extract_position(
              mutation
            )
          ) %>%
          
          filter(
            position == input$selected_position
          ) %>%
          
          arrange(
            desc(NBES)
          )
        
      })
      
      output$dataset_stage <- renderText({
        
        dataset_stage()
        
      })
      
      outputOptions(
        output,
        "dataset_stage",
        suspendWhenHidden = FALSE
      )
      
      # ----------------------------------------------------------------------
      # Dataset Information
      # ----------------------------------------------------------------------
      
      output$dataset_info <- renderPrint({
        
        ds <- current_dataset()
        
        cat("Name:", ds$name, "\n")
        cat("Stage:", ds$stage, "\n")
        cat("Source:", ds$source, "\n")
        cat("Created:", ds$created, "\n")
        
        if(ds$stage == "enrichment") {
          
          cat("\nEnrichment Metadata\n")
          cat("-------------------\n")
          
          cat(
            "Pre-sort Dataset:",
            ds$metadata$presort_dataset,
            "\n"
          )
          
          cat(
            "Selected Dataset:",
            ds$metadata$selected_dataset,
            "\n"
          )
          
        }
        
        if(ds$stage == "nbes") {
          
          cat("\nNBES Metadata\n")
          cat("-------------------\n")
          
          cat(
            "High Dataset:",
            ds$metadata$high_dataset,
            "\n"
          )
          
          cat(
            "Low Dataset:",
            ds$metadata$low_dataset,
            "\n"
          )
          
        }
        
        if (!is.null(ds$metadata)) {
          
          if (!is.null(ds$metadata$min_count)) {
            
            cat(
              "Min Count:",
              ds$metadata$min_count,
              "\n"
            )
            
          }
          
          if (!is.null(ds$metadata$min_freq)) {
            
            cat(
              "Min Frequency:",
              ds$metadata$min_freq,
              "\n"
            )
            
          }
          
        }
        
      })
      
      # ----------------------------------------------------------------------
      # Enrichment Summary
      # ----------------------------------------------------------------------
      
      output$enrichment_summary <- renderPrint({
        
        ds <- current_dataset()
        
        req(ds$stage == "enrichment")
        
        s <- ds$metadata$enrichment_summary
        
        cat(
          "Pre-sort variants:",
          s$n_presort_variants,
          "\n"
        )
        
        cat(
          "Selected variants:",
          s$n_selected_variants,
          "\n"
        )
        
        cat(
          "Variants used:",
          s$n_variants_used,
          "\n"
        )
        
        cat(
          "Variants dropped:",
          s$n_variants_dropped,
          "\n"
        )
        
        cat(
          "Percent dropped:",
          s$pct_variants_dropped,
          "%\n"
        )
        
        cat(
          "WT log2E:",
          round(
            s$wt_log2E,
            3
          ),
          "\n"
        )
        
        cat(
          "Pseudocount:",
          signif(
            s$pseudocount,
            3
          ),
          "\n"
        )
        
      })
      
      # ----------------------------------------------------------------------
      # NBES Summary
      # ----------------------------------------------------------------------
      
      output$nbes_summary <- renderPrint({
        
        ds <- current_dataset()
        
        req(ds$stage == "nbes")
        
        s <- ds$metadata$nbes_summary
        
        cat(
          "High Binder Dataset:",
          ds$metadata$high_dataset,
          "\n"
        )
        
        cat(
          "Low Binder Dataset:",
          ds$metadata$low_dataset,
          "\n\n"
        )
        
        cat(
          "High Variants:",
          s$n_high_variants,
          "\n"
        )
        
        cat(
          "Low Variants:",
          s$n_low_variants,
          "\n"
        )
        
        cat(
          "Variants Used:",
          s$n_variants_used,
          "\n"
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Mutation QC Metrics
      # ----------------------------------------------------------------------
      
      output$qc <- renderTable({
        
        req(
          current_dataset()$stage == "mutation"
        )
        
        as.data.frame(
          current_dataset()$results$qc
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Mutation QC Plots
      # ----------------------------------------------------------------------
      
      output$haplotype_distribution <- renderPlot({
        
        ds <- current_dataset()
        
        req(ds$stage == "mutation")
        
        min_count <- DEFAULT_MIN_COUNT
        
        if(
          !is.null(ds$metadata) &&
          !is.null(ds$metadata$min_count)
        ) {
          
          min_count <- ds$metadata$min_count
          
        }
        
        plot_haplotype_distribution(
          
          haplo_df_raw =
            ds$results$haplo_df_raw,
          
          min_count =
            min_count
          
        )
        
      })
      
      output$variant_classes <- renderPlot({
        
        ds <- current_dataset()
        
        req(ds$stage == "mutation")
        
        plot_variant_class_distribution(
          ds$results$haplo_df
        )
        
      })
      
      output$variant_reads <- renderPlot({
        
        ds <- current_dataset()
        
        req(ds$stage == "mutation")
        
        plot_variant_class_reads(
          ds$results$haplo_df
        )
        
      })
      
      # ----------------------------------------------------------------------
      # Mutation QC Table
      # ----------------------------------------------------------------------
      
      output$variant_qc <- DT::renderDT({
        
        ds <- current_dataset()
        
        req(ds$stage == "mutation")
        
        variant_qc_summary(
          ds$results$haplo_df
        )
        
      },
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # Enrichment Plots
      # ----------------------------------------------------------------------
      
      output$enrichment_distribution <- renderPlot({
        
        plot_enrichment_distribution(
          enrichment_df()
        )
        
      })
      
      output$position_enrichment <- renderPlot({
        
        plot_position_enrichment(
          position_summary()
        )
        
      })
      
      output$enrichment_heatmap <- renderPlot({
        
        ht <- plot_enrichment_heatmap(
          enrichment_df()
        )
        
        draw(ht)
        
      })
      
      # ----------------------------------------------------------------------
      # Enrichment Tables
      # ----------------------------------------------------------------------
      
      output$top_enriched <- DT::renderDT({
        
        top_enriched_variants(
          enrichment_df()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      output$top_depleted <- DT::renderDT({
        
        top_depleted_variants(
          enrichment_df()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      output$position_enrichment_summary <- DT::renderDT({
        
        position_summary_enrichment(
          current_dataset()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # NBES Plots
      # ----------------------------------------------------------------------
      
      output$nbes_distribution <- renderPlot({
        
        plot_nbes_distribution(
          nbes_df()
        )
        
      })
      
      output$position_nbes <- renderPlot({
        
        plot_position_nbes(
          position_summary()
        )
        
      })
      output$nbes_heatmap <- renderPlot({
        
        ds <- current_dataset()
        
        ht <- plot_nbes_heatmap(
          nbes_df(),
          ref_protein = ds$metadata$ref_protein
        )
        draw(ht)
        
      })
      
      output$position_variants <- DT::renderDT({
        
        position_variants()
        
      },
      options = list(
        pageLength = 25,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # NBES Tables
      # ----------------------------------------------------------------------
      
      output$top_positive_nbes <- DT::renderDT({
        
        top_positive_nbes(
          nbes_df()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      output$top_negative_nbes <- DT::renderDT({
        
        top_negative_nbes(
          nbes_df()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      output$position_nbes_summary <- DT::renderDT({
        
        position_summary_nbes(
          current_dataset()
        )
        
      },
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ))
      
      # ----------------------------------------------------------------------
      # Variant Browser
      # ----------------------------------------------------------------------
      
      variant_table <- reactive({
        
        ds <- current_dataset()
        
        if(ds$stage == "nbes") {
          
          return(
            ds$results$nbes
          )
          
        }
        
        if(ds$stage == "enrichment") {
          
          return(
            ds$results$enrichment
          )
          
        }
        
        switch(
          
          input$dataset_view,
          
          "All Variants" =
            ds$results$haplo_df,
          
          "Single Mutants" =
            ds$results$single_mutants,
          
          "Missense" =
            ds$results$missense,
          
          "Indels" =
            ds$results$indels,
          
          "Frameshifts" =
            ds$results$frameshifts,
          
          ds$results$single_mutants
          
        )
        
      })
      
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