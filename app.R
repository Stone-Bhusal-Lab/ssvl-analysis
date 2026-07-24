source("R/package_setup.R")

source("config/defaults.R")

source("R/datasets.R")
source("R/fastq_processing.R")
source("R/plots.R")
source("R/colours.R")
source("R/tables.R")
source("R/position_analysis.R")
source("R/enrichment_analysis.R")
source("R/nbes_analysis.R")

source("modules/mod_data_manager.R")
source("modules/mod_mutation_analysis.R")
source("modules/mod_visualisation.R")
source("modules/mod_comparison.R")
source("modules/mod_enrichment.R")
source("modules/mod_nbes.R")

check_required_packages()

load_required_packages()

# 8 GB file limit
options(shiny.maxRequestSize = 8 * 1024^3)


ui <- fluidPage(
  
  titlePanel("Stone Lab SSVL Workbench"),
  
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel(
      "Data Manager",
      mod_data_manager_ui("data")
    ),
    
    tabPanel(
      "Mutation Analysis",
      mod_mutation_analysis_ui("mutation")
    ),
    
    tabPanel(
      "Enrichment Analysis",
      mod_enrichment_ui("enrich")
    ),
    
    tabPanel(
      "NBES",
      mod_nbes_ui("nbes")
    ),
    
    tabPanel(
      "Visualisation",
      mod_visualisation_ui("viz")
    ),
    
    tabPanel(
      "Comparison",
      mod_comparison_ui("compare")
    )
  )
)

server <- function(input, output, session) {
  
  datasets <- reactiveVal(list())
  
  active_dataset <- reactiveVal(NULL)
  
  mod_data_manager_server(
    "data",
    datasets,
    active_dataset
  )
  
  mod_mutation_analysis_server(
    "mutation",
    datasets,
    active_dataset
  )
  
  mod_visualisation_server(
    "viz",
    datasets,
    active_dataset
  )
  
  mod_comparison_server(
    "compare",
    datasets,
    active_dataset
  )
  
  mod_enrichment_server(
    "enrich",
    datasets,
    active_dataset
  )
  
  mod_nbes_server(
    "nbes",
    datasets,
    active_dataset
  )
}

shinyApp(ui, server)