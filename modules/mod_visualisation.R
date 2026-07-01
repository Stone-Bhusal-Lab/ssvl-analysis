# UI
mod_visualisation_ui <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    h3("Visualisation"),
    
    DT::DTOutput(
      ns("variants")
    )
  )
}
# Server
mod_visualisation_server <- function(id,
                                     datasets,
                                     active_dataset){
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    current_dataset <- reactive({
      
      req(active_dataset())
      
      datasets()[[active_dataset()]]
    })
    
    output$variants <- DT::renderDT({
      
      req(current_dataset())
      
      current_dataset()$results$single_mutants
      
    })
  })
}
