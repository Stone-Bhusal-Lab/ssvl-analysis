library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(ggrepel)
library(stringr)
library(DT)
library(r3dmol)
library(htmlwidgets)

source("config/defaults.R")
source("R/fastq_processing.R")
source("R/plots.R")

app_password <- "test123"

taylor_muted <- c(
  A = "#99b300", R = "#3b3bb3", N = "#9933b3", D = "#cc3333",
  C = "#cccc33", Q = "#cc3399", E = "#cc3366", G = "#cc7a00",
  H = "#3366cc", I = "#66cc00", L = "#4dcc00", K = "#6633cc",
  M = "#33cc33", F = "#33cc66", P = "#cc9900", S = "#cc3300",
  T = "#cc6600", W = "#3399cc", Y = "#33cc99", V = "#80cc00"
)

colour_map <- c(taylor_muted, del = "black")

ui <- fluidPage(
  
  tags$head(
    tags$script(HTML("
      Shiny.addCustomMessageHandler('clicked_residue', function(x) {
        Shiny.setInputValue('pdb_click', x, {priority: 'event'});
      });
    "))
  ),
  
  conditionalPanel(
    condition = "!output.auth",
    fluidRow(column(4, offset = 4,
                    passwordInput("password", "Password"),
                    actionButton("login", "Login")
    ))
  ),
  
  conditionalPanel(
    condition = "output.auth",
    
    tabsetPanel(
      
      tabPanel(
        "FASTQ Analysis",
        
        br(),
        
        sidebarLayout(
          
          sidebarPanel(
            
            fileInput(
              "fastq",
              "Upload merged FASTQ"
            ),
            
            numericInput(
              "min_count",
              "Minimum read count",
              value = DEFAULT_MIN_COUNT
            ),
            
            numericInput(
              "min_freq",
              "Minimum frequency",
              value = DEFAULT_MIN_FREQ
            ),
            
            textAreaInput(
              "ref_protein",
              "Reference protein",
              value = DEFAULT_REF_PROTEIN,
              rows = 4
            ),
            
            textInput(
              "left_flank",
              "Left flank",
              value = DEFAULT_LEFT_FLANK
            ),
            
            textInput(
              "right_flank",
              "Right flank",
              value = DEFAULT_RIGHT_FLANK
            ),
            
            actionButton(
              "run_analysis",
              "Run Analysis"
            ),
            hr(),
            
            downloadButton(
              "download_single_mutants",
              "Download Single Mutants"
            ),
            
            downloadButton(
              "download_all_variants",
              "Download All Variants (Unfiltered)"
            )
            
          ),
          
          mainPanel(
            
            h3("Mutation Analysis"),
            fluidRow(
              column(3, wellPanel(
                h4("Reads"),
                textOutput("qc_total_reads")
              )),
              
              column(3, wellPanel(
                h4("Valid Inserts"),
                textOutput("qc_valid_inserts")
              )),
              
              column(3, wellPanel(
                h4("Extraction Rate"),
                textOutput("qc_extraction_rate")
              )),
              
              column(3, wellPanel(
                h4("Filtered Haplotypes"),
                textOutput("qc_filtered_haplotypes")
              ))
            ),
            plotOutput("haplotype_distribution"),
            DTOutput("analysis_table")
            
          )
        )
      ),
      
      tabPanel(
        "Explorer",
        
        titlePanel("Mutation Explorer + Structure"),
        
        sidebarLayout(
          
          sidebarPanel(
            width = 3,
            
            fileInput("files", "Upload datasets", multiple = TRUE),
            selectInput("dataset", "Dataset", choices = NULL),
            
            fileInput("pdb", "Upload PDB"),
            selectInput("chain", "Chain", choices = c("A","B","C")),
            
            sliderInput(
              "nbes_high_q",
              "NBES high",
              0.8,
              0.99,
              0.95
            ),
            
            sliderInput(
              "nbes_low_q",
              "NBES low",
              0.01,
              0.2,
              0.05
            ),
            
            selectInput(
              "structure_colour",
              "Structure colour",
              c("NBES","Deviation")
            ),
            
            actionButton(
              "clear_selection",
              "Clear selection"
            ),
            
            downloadButton(
              "plot",
              "Export plot"
            )
          ),
          
          mainPanel(
            
            fluidRow(
              column(6, plotlyOutput("interactive")),
              column(6, plotOutput("labelled"))
            ),
            
            r3dmolOutput(
              "structure",
              height = "500px"
            ),
            
            
            DTOutput("table")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  analysis_results <- eventReactive(
    input$run_analysis,
    {
      
      req(input$fastq)
      
      withProgress(
        message = "Processing FASTQ",
        value = 0,
        {
          
          incProgress(
            0.2,
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
            0.8,
            detail = "Finalising results"
          )
          
          result
        }
      )
    }
  )
  output$analysis_table <- renderDT({
    
    req(analysis_results())
    
    analysis_results()$single_mutants
    
  })
  
  output$haplotype_distribution <- renderPlot({
    
    req(analysis_results())
    
    plot_haplotype_distribution(
      haplo_df_raw = analysis_results()$haplo_df_raw,
      min_count = input$min_count
    )
    
  })
  output$qc_total_reads <- renderText({
    
    req(analysis_results())
    
    format(
      analysis_results()$qc$total_reads,
      big.mark = ","
    )
    
  })
  output$qc_valid_inserts <- renderText({
    
    req(analysis_results())
    
    format(
      analysis_results()$qc$valid_inserts,
      big.mark = ","
    )
    
  })
  output$qc_extraction_rate <- renderText({
    
    req(analysis_results())
    
    qc <- analysis_results()$qc
    
    paste0(
      round(
        100 * qc$valid_inserts / qc$total_reads,
        1
      ),
      "%"
    )
    
  })
  output$qc_filtered_haplotypes <- renderText({
    
    req(analysis_results())
    
    format(
      analysis_results()$qc$unique_haplotypes_filtered,
      big.mark = ","
    )
    
  })
  
  output$download_single_mutants <- downloadHandler(
    
    filename = function() {
      paste0(
        "single_mutants_",
        Sys.Date(),
        ".tsv"
      )
    },
    
    content = function(file) {
      
      req(analysis_results())
      
      write.table(
        analysis_results()$single_mutants,
        file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  output$download_all_variants <- downloadHandler(
    
    filename = function() {
      paste0(
        "all_variants_unfiltered_",
        Sys.Date(),
        ".tsv"
      )
    },
    
    content = function(file) {
      
      req(analysis_results())
      
      write.table(
        analysis_results()$haplo_df_raw,
        file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  
  output$analysis_summary <- renderPrint({
    
    req(analysis_results())
    
    qc <- analysis_results()$qc
    
    cat(
      "Total reads:",
      qc$total_reads,
      "\n"
    )
    
    cat(
      "Valid inserts:",
      qc$valid_inserts,
      "\n"
    )
    
    cat(
      "Extraction rate:",
      round(
        100 * qc$valid_inserts / qc$total_reads,
        1
      ),
      "%\n"
    )
    
    cat(
      "Unique haplotypes (raw):",
      qc$unique_haplotypes_raw,
      "\n"
    )
    
    cat(
      "Unique haplotypes (filtered):",
      qc$unique_haplotypes_filtered,
      "\n\n"
    )
    
    cat(
      "Single mutant entries:",
      nrow(analysis_results()$single_mutants),
      "\n"
    )
    
    cat(
      "Missense variants:",
      nrow(analysis_results()$missense),
      "\n"
    )
    
    cat(
      "In-frame deletions:",
      nrow(analysis_results()$indels),
      "\n"
    )
    
    cat(
      "Frameshifts:",
      nrow(analysis_results()$frameshifts),
      "\n"
    )
    
  })
  
  
  auth <- reactiveVal(FALSE)
  observeEvent(input$login, {
    if (input$password == app_password) {
      auth(TRUE)
    }
  })
  output$auth <- reactive(auth())
  outputOptions(output, "auth", suspendWhenHidden = FALSE)
  
  datasets <- reactiveVal(list())
  
  observeEvent(input$files,{
    new <- lapply(input$files$datapath, read.csv)
    names(new) <- input$files$name
    datasets(c(datasets(), new))
  })
  
  observe({
    updateSelectInput(session, "dataset", choices = names(datasets()))
  })
  
  df_proc <- reactive({
    req(input$dataset)
    datasets()[[input$dataset]] %>%
      mutate(
        position = as.numeric(str_extract(mutation, "\\d+")),
        alt = str_extract(mutation, "[A-Z]+$")
      )
  })
  
  selected <- reactiveVal(character(0))
  
  classified <- reactive({ df_proc() })
  
  # ---------------- INTERACTIVE PLOT
  output$interactive <- renderPlotly({
    
    dfc <- classified()
    
    p <- ggplot(dfc,
                aes(NBES, position_deviation, key = mutation,
                    text = paste("Residue:", position))
    ) +
      geom_point(aes(colour = alt)) +
      scale_colour_manual(values = colour_map)
    
    ggplotly(p, source = "sel", tooltip = "text")
  })
  
  # Plot selections
  observeEvent(event_data("plotly_click", source="sel"), {
    
    ed <- event_data("plotly_click", source="sel")
    if (is.null(ed)) return()
    
    cur <- selected()
    if (ed$key %in% cur) {
      selected(setdiff(cur, ed$key))
    } else {
      selected(c(cur, ed$key))
    }
  })
  
  observeEvent(input$clear_selection,{
    selected(character(0))
  })
  
  # ---------------- LABELLED
  labelled_plot <- reactive({
    
    dfc <- classified()
    
    p <- ggplot(dfc, aes(NBES, position_deviation)) +
      geom_point(aes(colour = alt)) +
      scale_colour_manual(values = colour_map)
    
    sel <- dfc %>% filter(mutation %in% selected())
    
    if (nrow(sel)>0) {
      p <- p +
        geom_text_repel(
          data = sel,
          aes(label = mutation),
          max.overlaps = Inf
        )
    }
    
    p
  })
  
  output$labelled <- renderPlot(labelled_plot())
  
  # ---------------- 3D VIEW WITH CLICK HOOK
  output$structure <- renderR3dmol({
    
    req(input$pdb)
    
    pdb_text <- paste(readLines(input$pdb$datapath), collapse="\n")
    
    r3dmol() %>%
      m_add_model(pdb_text, "pdb") %>%
      m_set_style(style = m_style_cartoon()) %>%
      m_add_js("
        function(el,x){
          var viewer = this;
          viewer.setClickable({}, true, function(atom){
            if(atom.resi){
              Shiny.setInputValue('pdb_click', atom.resi);
            }
          });
        }
      ")
  })
  
  # ---------------- STRUCTURE → PLOT
  observeEvent(input$pdb_click, {
    
    req(input$pdb_click)
    
    dfc <- classified()
    
    muts <- dfc %>%
      filter(position == input$pdb_click) %>%
      pull(mutation)
    
    selected(unique(c(selected(), muts)))
  })
  
  # ---------------- PLOT → STRUCTURE
  observe({
    
    req(input$pdb)
    
    dfc <- classified()
    proxy <- r3dmol_proxy("structure")
    
    proxy %>% m_set_style(style = m_style_cartoon(color="grey"))
    
    sel <- dfc %>%
      filter(mutation %in% selected()) %>%
      pull(position)
    
    if(length(sel)>0){
      proxy %>%
        m_set_style(
          sel = list(resi = sel),
          style = m_style_stick(color="yellow")
        )
    }
  })
  
  # ---------------- EXPORT
  output$plot <- downloadHandler(
    filename=function() "plot.png",
    content=function(file){
      ggsave(file, labelled_plot(), width=8,height=6)
    }
  )
  
  output$table <- renderDT(classified())
}

shinyApp(ui, server)