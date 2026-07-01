# UI
mod_data_manager_ui <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    fluidRow(
      
      column(
        4,
        
        h3("Load Dataset"),
        
        fileInput(
          ns("load_rds"),
          "Load RDS Dataset"
        )
      ),
      
      column(
        8,
        
        h3("Dataset Library"),
        
        tableOutput(
          ns("dataset_table")
        )
      )
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        6,
        
        selectInput(
          ns("active_dataset"),
          "Active Dataset",
          choices = NULL
        )
      ),
      
      column(
        6,
        
        downloadButton(
          ns("save_dataset"),
          "Export Active Dataset"
        )
      )
    )
  )
}
# Server
mod_data_manager_server <- function(id,
                                    datasets,
                                    active_dataset){
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    observeEvent(input$load_rds, {
      
      ds <- readRDS(
        input$load_rds$datapath
      )
      
      current <- datasets()
      
      current[[ds$name]] <- ds
      
      datasets(current)
      
    })
    
    observe({
      
      updateSelectInput(
        session,
        "active_dataset",
        choices = names(datasets())
      )
    })
    
    observeEvent(
      input$active_dataset,
      active_dataset(input$active_dataset)
    )
    
    output$dataset_table <- renderTable({
      
      req(length(datasets()) > 0)
      
      do.call(
        rbind,
        lapply(
          datasets(),
          function(x){
            data.frame(
              Name = x$name,
              Stage = x$stage,
              Source = x$source
            )
          }
        )
      )
      
    })
    
    output$save_dataset <- downloadHandler(
      
      filename = function(){
        
        paste0(
          active_dataset(),
          ".rds"
        )
      },
      
      content = function(file){
        
        saveRDS(
          datasets()[[active_dataset()]],
          file
        )
      }
    )
  })
}
