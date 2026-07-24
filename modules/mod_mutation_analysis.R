# UI
mod_mutation_analysis_ui <- function(id){
  
  ns <- NS(id)
  
  sidebarLayout(
    
    sidebarPanel(
      
      textInput(
        ns("dataset_name"),
        "Dataset Name"
      ),
      
      fileInput(
        ns("fastq"),
        "FASTQ"
      ),
      
      helpText(
        "Large FASTQ files may take several minutes to upload and process."
      ),
      
      numericInput(
        ns("min_count"),
        "Minimum Count",
        DEFAULT_MIN_COUNT
      ),
      
      numericInput(
        ns("min_freq"),
        "Minimum Frequency",
        DEFAULT_MIN_FREQ
      ),
      
      textAreaInput(
        ns("ref_protein"),
        "Reference Protein",
        DEFAULT_REF_PROTEIN
      ),
      
      textInput(
        ns("left_flank"),
        "Left Flank",
        DEFAULT_LEFT_FLANK
      ),
      
      textInput(
        ns("right_flank"),
        "Right Flank",
        DEFAULT_RIGHT_FLANK
      ),
      
      actionButton(
        ns("run"),
        "Run Analysis"
      )
    ),
    
    mainPanel(
      
      tableOutput(
        ns("qc")
      ),
      
      fluidRow(
        
        column(
          4,
          plotOutput(
            ns("distribution")
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
      )
      
    )
  )
}

# Server
mod_mutation_analysis_server <- function(id,
                                         datasets,
                                         active_dataset){
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    results <- eventReactive(
      input$run,
      {
        
        req(input$fastq)
        
        withProgress(
          
          message = "Processing FASTQ",
          value = 0,
          
          {
            
            incProgress(
              0.1,
              detail = "Reading sequences"
            )
            
            result <- run_mutation_analysis(
              fastq_file = input$fastq$datapath,
              ref_protein = input$ref_protein,
              left_flank = input$left_flank,
              right_flank = input$right_flank,
              min_count = input$min_count,
              min_freq = input$min_freq
            )
            
            incProgress(
              0.9,
              detail = "Finalising results"
            )
            
            result
            
          }
          
        )
        
      }
    )
    
    observeEvent(results(), {
      
      ds <- list(
        
        name = input$dataset_name,
        
        stage = "mutation",
        
        source = "fastq",
        
        created = Sys.time(),
        
        metadata = list(
          
          ref_protein = input$ref_protein,
          
          left_flank = input$left_flank,
          
          right_flank = input$right_flank,
          
          min_count = input$min_count,
          
          min_freq = input$min_freq
          
        ),
        
        results = results()
        
      )
      
      current <- datasets()
      
      current[[ds$name]] <- ds
      
      datasets(current)
      
      active_dataset(ds$name)
      
      showNotification(
        paste(
          "Dataset created:",
          ds$name
        ),
        type = "message",
        duration = 5
      )
      
      updateTabsetPanel(
        session,
        "main_tabs",
        selected = "Visualisation"
      )
      
    })
    
    output$qc <- renderTable({
      
      req(results())
      
      as.data.frame(
        results()$qc
      )
      
    })
    
    output$distribution <- renderPlot({
      
      req(results())
      
      plot_haplotype_distribution(
        haplo_df_raw = results()$haplo_df_raw,
        min_count = input$min_count
      )
      
    })
    
    output$variant_classes <- renderPlot({
      
      req(results())
      
      plot_variant_class_distribution(
        results()$haplo_df_raw
      )
      
    })
    
    output$variant_reads <- renderPlot({
      
      req(results())
      
      plot_variant_class_reads(
        results()$haplo_df_raw
      )
      
    })
    
    output$variant_qc <- DT::renderDT({
      
      req(results())
      
      variant_qc_summary(
        results()$haplo_df_raw
      )
      
    },
    options = list(
      pageLength = 10,
      scrollX = TRUE
    ))
    
  })
}