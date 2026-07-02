# ============================================================================
# DATA MANAGER MODULE - UI
# ============================================================================

mod_data_manager_ui <- function(id) {
  
  ns <- NS(id)
  
  fluidPage(
    
    # ------------------------------------------------------------------------
    # Dataset Import & Library
    # ------------------------------------------------------------------------
    
    fluidRow(
      
      column(
        
        width = 4,
        
        h3("Dataset Management"),
        
        fileInput(
          ns("load_rds"),
          "Load Dataset(s) (.rds)",
          accept = ".rds",
          multiple = TRUE
        ),
        
        br(),
        
        downloadButton(
          ns("save_dataset"),
          "Export Active Dataset"
        ),
        
        hr(),
        
        selectInput(
          ns("active_dataset"),
          "Active Dataset",
          choices = NULL
        ),
        
        hr(),
        
        textInput(
          ns("new_name"),
          "Rename Active Dataset"
        ),
        
        actionButton(
          ns("rename_dataset"),
          "Rename Dataset"
        ),
        
        br(),
        br(),
        
        actionButton(
          ns("duplicate_dataset"),
          "Duplicate Dataset"
        ),
        
        br(),
        br(),
        
        actionButton(
          ns("delete_dataset"),
          "Delete Dataset"
        )
        
      ),
      
      column(
        
        width = 8,
        
        h3("Dataset Library"),
        
        DT::DTOutput(
          ns("dataset_table")
        )
        
      )
      
    ),
    
    hr(),
    
    # ------------------------------------------------------------------------
    # Active Dataset Summary
    # ------------------------------------------------------------------------
    
    h3("Dataset Summary"),
    
    verbatimTextOutput(
      ns("active_dataset_display")
    )
    
  )
  
}


# ============================================================================
# DATA MANAGER MODULE - SERVER
# ============================================================================

mod_data_manager_server <- function(
    id,
    datasets,
    active_dataset
) {
  
  moduleServer(
    
    id,
    
    function(input, output, session) {
      
      # ----------------------------------------------------------------------
      # Load RDS Dataset
      # ----------------------------------------------------------------------
      
      observeEvent(input$load_rds, {
        
        req(input$load_rds)
        
        current <- datasets()
        
        loaded_names <- character()
        
        for (i in seq_along(input$load_rds$datapath)) {
          
          ds <- readRDS(
            input$load_rds$datapath[i]
          )
          
          name <- ds$name
          
          # Handle duplicate names
          
          if (name %in% names(current)) {
            
            counter <- 1
            
            candidate <- paste0(name, "_", counter)
            
            while (candidate %in% names(current)) {
              
              counter <- counter + 1
              
              candidate <- paste0(
                name,
                "_",
                counter
              )
              
            }
            
            name <- candidate
            
            ds$name <- name
            
          }
          
          current[[name]] <- ds
          
          loaded_names <- c(
            loaded_names,
            name
          )
          
        }
        
        datasets(current)
        
        if (length(loaded_names) > 0) {
          
          active_dataset(
            loaded_names[1]
          )
          
        }
        
        showNotification(
          paste(
            "Loaded",
            length(loaded_names),
            "dataset(s)"
          ),
          type = "message"
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Update Dataset Selector
      # ----------------------------------------------------------------------
      
      observe({
        
        updateSelectInput(
          
          session,
          
          "active_dataset",
          
          choices = names(datasets()),
          
          selected = active_dataset()
          
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Manual Dataset Selection
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$active_dataset,
        {
          
          active_dataset(
            input$active_dataset
          )
          
        }
      )
      
      
      
      # ----------------------------------------------------------------------
      # Dataset Summary Table
      # ----------------------------------------------------------------------
      
      dataset_summary <- reactive({
        
        if (length(datasets()) == 0) {
          
          return(
            data.frame(
              Status = "No datasets loaded"
            )
          )
          
        }
        
        do.call(
          
          rbind,
          
          lapply(
            
            datasets(),
            
            function(ds) {
              
              if (ds$stage == "mutation") {
                
                qc <- ds$results$qc
                
                return(
                  
                  data.frame(
                    
                    Dataset = ds$name,
                    
                    Stage = ds$stage,
                    
                    Source = ds$source,
                    
                    Reads = qc$total_reads,
                    
                    Variants = nrow(
                      ds$results$single_mutants
                    ),
                    
                    Created = format(
                      ds$created,
                      "%Y-%m-%d %H:%M"
                    ),
                    
                    check.names = FALSE
                    
                  )
                  
                )
                
              }
              
              if (ds$stage == "enrichment") {
                
                s <- ds$metadata$enrichment_summary
                
                return(
                  
                  data.frame(
                    
                    Dataset = ds$name,
                    
                    Stage = ds$stage,
                    
                    Source = ds$source,
                    
                    Reads = NA,
                    
                    Variants = s$n_variants_used,
                    
                    Created = format(
                      ds$created,
                      "%Y-%m-%d %H:%M"
                    ),
                    
                    check.names = FALSE
                    
                  )
                  
                )
                
              }
              
              if (ds$stage == "nbes") {
                
                s <- ds$metadata$nbes_summary
                
                return(
                  
                  data.frame(
                    
                    Dataset = ds$name,
                    
                    Stage = ds$stage,
                    
                    Source = ds$source,
                    
                    Reads = NA,
                    
                    Variants = s$n_variants_used,
                    
                    Created = format(
                      ds$created,
                      "%Y-%m-%d %H:%M"
                    ),
                    
                    check.names = FALSE
                    
                  )
                  
                )
                
              }
              
            }
            
          )
          
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Render Dataset Table
      # ----------------------------------------------------------------------
      
      output$dataset_table <- DT::renderDT({
        
        DT::datatable(
          
          dataset_summary(),
          
          selection = "single",
          
          options = list(
            pageLength = 10,
            scrollX = TRUE
          ),
          
          rownames = FALSE
          
        )
        
      })
      
      
      
      # ----------------------------------------------------------------------
      # Activate Dataset by Clicking Table Row
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$dataset_table_rows_selected,
        {
          
          req(
            input$dataset_table_rows_selected
          )
          
          idx <- input$dataset_table_rows_selected
          
          selected_name <-
            dataset_summary()$Dataset[idx]
          
          active_dataset(
            selected_name
          )
          
        }
      )
      
      
      
      # ----------------------------------------------------------------------
      # Rename Dataset
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$rename_dataset,
        {
          
          req(active_dataset())
          req(input$new_name)
          
          old_name <- active_dataset()
          new_name <- trimws(input$new_name)
          
          req(nchar(new_name) > 0)
          
          current <- datasets()
          
          current[[new_name]] <-
            current[[old_name]]
          
          current[[new_name]]$name <-
            new_name
          
          current[[old_name]] <- NULL
          
          datasets(current)
          
          active_dataset(new_name)
          
          showNotification(
            paste(
              "Renamed dataset to",
              new_name
            ),
            type = "message"
          )
          
        }
      )
      
      
      
      # ----------------------------------------------------------------------
      # Duplicate Dataset
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$duplicate_dataset,
        {
          
          req(active_dataset())
          
          current <- datasets()
          
          ds <- current[[active_dataset()]]
          
          new_name <- paste0(
            ds$name,
            "_copy"
          )
          
          counter <- 1
          
          while (
            new_name %in% names(current)
          ) {
            
            new_name <- paste0(
              ds$name,
              "_copy_",
              counter
            )
            
            counter <- counter + 1
            
          }
          
          ds$name <- new_name
          
          current[[new_name]] <- ds
          
          datasets(current)
          
          showNotification(
            paste(
              "Created",
              new_name
            ),
            type = "message"
          )
          
        }
      )
      
      
      
      # ----------------------------------------------------------------------
      # Delete Dataset
      # ----------------------------------------------------------------------
      
      observeEvent(
        input$delete_dataset,
        {
          
          req(active_dataset())
          
          deleted_name <- active_dataset()
          
          current <- datasets()
          
          current[[deleted_name]] <- NULL
          
          datasets(current)
          
          if (length(current) > 0) {
            
            active_dataset(
              names(current)[1]
            )
            
          } else {
            
            active_dataset(NULL)
            
          }
          
          showNotification(
            paste(
              "Deleted",
              deleted_name
            ),
            type = "warning"
          )
          
        }
      )
      
      
      
      # ----------------------------------------------------------------------
      # Active Dataset Summary
      # ----------------------------------------------------------------------
      
      output$active_dataset_display <- renderPrint({
        
        req(active_dataset())
        
        ds <- datasets()[[active_dataset()]]
        
        cat(
          "Dataset:",
          ds$name,
          "\n"
        )
        
        cat(
          "Stage:",
          ds$stage,
          "\n"
        )
        
        cat(
          "Source:",
          ds$source,
          "\n"
        )
        
        cat(
          "Created:",
          ds$created,
          "\n\n"
        )
        
        if(ds$stage == "mutation"){
          
          qc <- ds$results$qc
          
          cat(
            "QC Summary\n"
          )
          
          cat(
            "----------\n"
          )
          
          cat(
            "Reads:",
            qc$total_reads,
            "\n"
          )
          
          cat(
            "Valid Inserts:",
            qc$valid_inserts,
            "\n"
          )
          
          cat(
            "Filtered Haplotypes:",
            qc$unique_haplotypes_filtered,
            "\n"
          )
          
        }
        
        if(ds$stage == "enrichment"){
          
          s <- ds$metadata$enrichment_summary
          
          cat(
            "Enrichment Summary\n"
          )
          
          cat(
            "------------------\n"
          )
          
          cat(
            "Variants Used:",
            s$n_variants_used,
            "\n"
          )
          
          cat(
            "Variants Dropped:",
            s$n_variants_dropped,
            "\n"
          )
          
          cat(
            "WT log2E:",
            round(
              s$wt_log2E,
              3
            ),
            "\n"
          )
          
        }
        
        if(ds$stage == "nbes"){
          
          s <- ds$metadata$nbes_summary
          
          cat(
            "NBES Summary\n"
          )
          
          cat(
            "------------\n"
          )
          
          cat(
            "Variants Used:",
            s$n_variants_used,
            "\n"
          )
          
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
        
      })
      # ----------------------------------------------------------------------
      # Export Active Dataset
      # ----------------------------------------------------------------------
      
      output$save_dataset <- downloadHandler(
        
        filename = function() {
          
          req(active_dataset())
          
          paste0(
            active_dataset(),
            ".rds"
          )
          
        },
        
        content = function(file) {
          
          req(active_dataset())
          
          saveRDS(
            
            datasets()[[active_dataset()]],
            
            file
            
          )
          
        }
        
      )
      
    }
    
  )
  
}