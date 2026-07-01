# UI
mod_enrichment_ui <- function(id){
  
  fluidPage(
    
    h3("Enrichment Analysis"),
    
    p(
      "Enrichment functionality coming next."
    )
  )
}

mod_enrichment_server <- function(id,
                                  datasets,
                                  active_dataset){
  
  moduleServer(id, function(input,
                            output,
                            session){
    
  })
}
# Server
