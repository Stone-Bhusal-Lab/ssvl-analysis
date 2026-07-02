# ============================================================================
# NBES MODULE - UI
# ============================================================================

mod_nbes_ui <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    sidebarLayout(
      
      sidebarPanel(
        
        h3("NBES Analysis"),
        
        selectInput(
          ns("high_dataset"),
          "High Binder Enrichment Dataset",
          choices = NULL
        ),
        
        selectInput(
          ns("low_dataset"),
          "Low Binder Enrichment Dataset",
          choices = NULL
        ),
        
        textInput(
          ns("dataset_name"),
          "NBES Dataset Name"
        ),
        
        actionButton(
          ns("run_nbes"),
          "Generate NBES"
        )
        
      ),
      
      mainPanel(
        
        h3("Dataset Preview"),
        
        tableOutput(
          ns("dataset_summary")
        ),
        
        hr(),
        
        DT::DTOutput(
          ns("preview")
        )
        
      )
      
    )
    
  )
  
}

# ============================================================================
# NBES MODULE - SERVER
# ============================================================================

mod_nbes_server <- function(
    id,
    datasets,
    active_dataset
){
  
  moduleServer(
    
    id,
    
    function(input, output, session){
      
      # ----------------------------------------------------------------------
      # Available Enrichment Datasets
      # ----------------------------------------------------------------------
      
      observe({
        
        enrichment_datasets <- names(
          
          Filter(
            
            function(x){
              
              !is.null(x$stage) &&
                x$stage == "enrichment"
              
            },
            
            datasets()
            
          )
          
        )
        
        updateSelectInput(
          session,
          "high_dataset",
          choices = enrichment_datasets
        )
        
        updateSelectInput(
          session,
          "low_dataset",
          choices = enrichment_datasets
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Dataset Reactives
      # ----------------------------------------------------------------------
      
      high_dataset <- reactive({
        
        req(input$high_dataset)
        
        datasets()[[input$high_dataset]]
        
      })
      
      
      
      low_dataset <- reactive({
        
        req(input$low_dataset)
        
        datasets()[[input$low_dataset]]
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Summary
      # ----------------------------------------------------------------------
      
      output$dataset_summary <- renderTable({
        
        req(
          high_dataset(),
          low_dataset()
        )
        
        data.frame(
          
          Dataset = c(
            high_dataset()$name,
            low_dataset()$name
          ),
          
          Role = c(
            "High Binder",
            "Low Binder"
          ),
          
          Variants = c(
            nrow(
              high_dataset()$results$enrichment
            ),
            nrow(
              low_dataset()$results$enrichment
            )
          )
          
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Preview
      # ----------------------------------------------------------------------
      
      output$preview <- DT::renderDT({
        
        req(high_dataset())
        
        high_dataset()$results$enrichment
        
      },
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ))
      
      
      
      # ----------------------------------------------------------------------
      # Run NBES Analysis
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$run_nbes,
        {
          
          validate(
            
            need(
              input$high_dataset != input$low_dataset,
              "High and low binder datasets must be different."
            )
            
          )
          
          
          
          nbes_results <- calculate_nbes(
            
            high_df =
              high_dataset()$results$enrichment,
            
            low_df =
              low_dataset()$results$enrichment
            
          )
          
          
          
          dataset_name <- trimws(
            input$dataset_name
          )
          
          if(dataset_name == ""){
            
            dataset_name <- paste0(
              
              high_dataset()$name,
              
              "_vs_",
              
              low_dataset()$name,
              
              "_NBES"
              
            )
            
          }
          
          
          
          nbes_ds <- create_nbes_dataset(
            
            nbes_results =
              nbes_results,
            
            dataset_name =
              dataset_name,
            
            high_dataset =
              input$high_dataset,
            
            low_dataset =
              input$low_dataset
            
          )
          
          
          
          current <- datasets()
          
          current[[dataset_name]] <- nbes_ds
          
          datasets(current)
          
          active_dataset(
            dataset_name
          )
          
          
          
          showNotification(
            
            paste(
              "Created NBES dataset:",
              dataset_name
            ),
            
            type = "message",
            duration = 5
            
          )
          
        }
        
      )
      
    }
    
  )
  
}