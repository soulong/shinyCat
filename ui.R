

file_accept_type <- c(".csv", ".xlsx")

meta_data_fixed_cols <- 1:4

jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"


##### ui #####
ui <- dashboardPage(
  #  skin = "black",
  
  ##### **hearder #####
  dashboardHeader(
    title = "shinyCat"
  ),
  
  ##### **sidebar #####
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload data", tabName="upload", selected = TRUE)
    ),
    
    sidebarMenu(
      menuItem("Sample QC", tabName="sample_qc")
    ), 
    
    sidebarMenu(
      menuItem("DGE analysis", tabName="dge_analysis")
    ), 
    
    sidebarMenu(
      menuItem("Enrichment",
               startExpanded = TRUE,
               menuSubItem("ORA", tabName = "ora"),
               menuSubItem("GSEA", tabName = "gsea")
      )
    )
    
  ), # end of dashboardSidebar
  
  ##### **body #####
  dashboardBody(
    useShinyjs(),
    #  extendShinyjs(text = jscode),
    
    tabItems(
      
      tabItem(tabName = "upload",
              fluidRow(
                box(id = "upload_box", title = "Upload data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    fileInput("counts_data",  paste0("Gene expression raw counts (", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
                    uiOutput("counts_data_select_sheet"),
                    fileInput("meta_data", paste0("Sample info metadata (optional, ", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
                    uiOutput("meta_data_select_sheet"),
                    #   actionButton("upload_data_submit", "Submit"),
                    disable(actionButton("upload_data_submit", "Submit"))
                    #   div(align="center", disable(actionButton("upload_data_submit", "Submit")))  
                ),
                # view counts
                box(title = "You can edit meta data here", status = "warning", solidHeader = TRUE, width = 6,
                    rHandsontableOutput("meta_data_table")
                ),
                # edit sample info
                box(title = "DESeq2 colData", status = "primary", solidHeader = TRUE, width = 6,
                    fluidRow(
                      column(9,
                             div(align = "center", uiOutput("choose_factors_ui"))),
                      column(3,
                             div(align = "center", actionButton("upload_data_apply", "Apply")))
                    ),
                    rHandsontableOutput("col_data_table")
                )
              )
      ),
      
      tabItem(tabName = "sample_qc",
      ),
      
      tabItem(tabName = "dge_analysis",
      ),
      
      tabItem(tabName = "ora",
      ),
      
      tabItem(tabName = "gsea",
              fluidRow(
                
              )
      )
      
    ) # end of tabItems
  ) # end of dashboardBody
  
) # end of dashboardPage

