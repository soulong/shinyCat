
##### UTILITIES #####
source("utilities.R")

##### MODULES #####
source("prepare_data.R")
source("sample_qc.R")
source("deseq2_compare.R")
source("enrichment_ora.R")


##### UI #####
ui <- dashboardPage(
  #  skin = "black",
  
  ##### **header #####
  dashboardHeader(
    title = "shinyCat"
  ),
  
  ##### **sidebar #####
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName="prepare_data", selected = TRUE)
    ),
    
    sidebarMenu(
      menuItem("Exploratory", tabName="sample_qc")
    ), 
    
    sidebarMenu(
      menuItem("DGE analysis", tabName = "deseq2_compare")
    ), 
    
    sidebarMenu(
      menuItem("Enrichment",
               startExpanded = TRUE,
               menuSubItem("ORA", tabName = "enrichment_ora"),
               menuSubItem("GSEA", tabName = "enrichment_gsea")
      )
    )
    
  ), # end of sidebar
  
  ##### **body #####
  dashboardBody(
    useShinyjs(),
    # extendShinyjs(text = jscode_collapse),
    
    tabItems(
      tabItem(tabName = "prepare_data",
              prepare_data_ui("prepare_data")
      ),
      tabItem(tabName = "sample_qc",
              sample_qc_ui("sample_qc")
      ),
      tabItem(tabName = "deseq2_compare",
              deseq2_compare_ui("deseq2_compare")
      ),
      tabItem(tabName = "enrichment_ora",
              enrichment_ora_ui("enrichment_ora")
      )
      # tabItem(tabName = "enrichment_gsea",
      #         enrichment_ora_ui("enrichment_gsea")
      # )
    )
  ) # end of body
)




##### SERVER #####
server <- function(input, output, session) {

  ##### upload #####
  data <- prepare_data_server("prepare_data")
  # observe({ print(data$coldata()) })
  
  ##### dds qc #####
  dds <- sample_qc_server("sample_qc", data$counts, data$coldata)
  
  ##### get results #####
  dds <- reactive({ read_rds("dds.rds") })
  observe({ glimpse(dds) })
  results <- deseq2_compare_server("deseq2_compare", dds)
  
  ##### enrichment #####
  results <- reactive({ read_rds("results.rds") })
  observe({ glimpse(results$results_sig()) })
  enrichment <- enrichment_ora_server("enrichment_ora", results$results_sig)
  
  
  
  # observe({ 
  #   validate(need(results()))
  #   glimpse(results()) })
  
  # observe({ glimpse(results()) })
  
  # # enable run_deseq2_apply buttom
  # observe({
  #   req(upload_counts(), col_data())
  #   enable("run_deseq2_apply")
  # })

}



##### RUN ##### 
shinyApp(ui, server)
