# UI
mod_comparison_ui <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    selectizeInput(
      ns("compare"),
      "Datasets",
      choices = NULL,
      multiple = TRUE
    ),
    
    tableOutput(
      ns("summary")
    )
  )
}
# Server
mod_comparison_server <- function(id,
                                  datasets,
                                  active_dataset){
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    observe({
      
      updateSelectizeInput(
        session,
        "compare",
        choices = names(datasets())
      )
    })
    
    output$summary <- renderTable({
      
      req(length(input$compare) > 0)
      
      do.call(
        rbind,
        lapply(
          input$compare,
          function(x){
            
            ds <- datasets()[[x]]
            
            data.frame(
              Dataset = x,
              Reads = ds$results$qc$total_reads,
              Valid = ds$results$qc$valid_inserts
            )
          }
        )
      )
    })
  })
}