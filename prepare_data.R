
# module_ui <- function(id, xyz) {
#   ns <- NS(id)
#   tagList(
#     fileInput(ns("id"), zyx),
#     ...
#   )
# }
prepare_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # upload file
      box(id = ns("upload_box"), title = "Upload data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
          fileInput(ns("counts_file"),  paste0("Gene expression raw counts (", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
          uiOutput(ns("counts_sheet_ui")),
          fileInput(ns("metadata_file"), paste0("Sample info metadata (optional, ", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
          uiOutput(ns("metadata_sheet_ui")),
          disabled(actionButton(ns("upload"), "Upload"))
      ),
      # view metadata
      box(title = "Edit meta data if necessary", status = "warning", solidHeader = TRUE, width = 6,
          rHandsontableOutput(ns("metadata_table"))
      ),
      # edit sample info
      box(title = "DESeq2 colData", status = "primary", solidHeader = TRUE, width = 6,
          # div(align = "center", uiOutput(ns("choose_factors_ui")) )
          fluidRow(
            column(9, div(align = "center", uiOutput(ns("choose_factors_ui")) )),
            # column(3,
            #        div(align = "center", disabled(actionButton(ns("deseq2_apply"), "Apply"))) )),
            # ),
            rHandsontableOutput(ns("coldata_table")) 
          )
      )
    )
  )
}



# module_server <- function(id, xyz, ...) {
#   moduleServer(
#     id,
#     function(input, output, session) {
#       fun1(xyz())
#       ...
#     })
# }
prepare_data_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ##### ** upload counts file #####
      counts_file <- reactive({
        validate(need(input$counts_file, message = "upload counts file"))
        input$counts_file })
      # file type
      counts_filetype <- reactive({
        if_else(file_extension(counts_file()$name) == "xlsx", "xlsx", "csv") })
      # select count data sheet UI
      observe({
        if(counts_filetype() == "xlsx") {
          counts_sheets <- excel_sheets(counts_file()$datapath)
          output$counts_sheet_ui <- renderUI({
            ns <- session$ns
            radioButtons(ns("counts_data_sheet_selected"), label = NULL,
                         choiceNames = counts_sheets,
                         choiceValues = counts_sheets,
                         selected = last(counts_sheets), inline = TRUE) })
        } else output$counts_sheet_ui <- NULL
      })
      # read counts data
      counts <- eventReactive(input$upload, {
        if(is.null(counts_file()$datapath)) return(NULL)
        showNotification("Read counts data ...", duration = 2, type = "message")
        if(counts_filetype() == "xlsx") {
          read_xlsx(counts_file()$datapath, sheet = input$counts_data_sheet_selected)
        } else read_csv(counts_file()$datapath)
      })
      
      
      
      ##### ** upload metadata file #####
      metadata_file <- reactive({
        validate(need(input$metadata_file, message = "upload metadata file"))
        input$metadata_file })
      # file type
      metadata_filetype <- reactive({
        if_else(file_extension(metadata_file()$name) == "xlsx", "xlsx", "csv") })
      # select metadata sheet UI
      observe({
        if(metadata_filetype() == "xlsx") {
          metadata_sheets <- excel_sheets(metadata_file()$datapath)
          output$metadata_sheet_ui <- renderUI({
            ns <- session$ns
            radioButtons(ns("metadata_sheet_selected"), label = NULL,
                         choiceNames = metadata_sheets,
                         choiceValues = metadata_sheets,
                         selected = first(metadata_sheets), inline = TRUE) })
        } else output$metadata_sheet_ui <- NULL
      })
      # read metadata
      metadata <- eventReactive(input$upload, {
        if(is.null(metadata_file()$datapath)) return(NULL)
        showNotification("Read metadata ...", duration = 1, type = "message")
        if(metadata_filetype() == "xlsx") {
          read_xlsx(metadata_file()$datapath, sheet = input$metadata_sheet_selected, col_types = "text")
        } else read_csv(metadata_file()$datapath, col_types = "text")
      })
      
      # enable button
      observe({
        req(counts_file(), metadata_file())
        enable("upload") })
      
      # collapse upload box
      # observeEvent(input$upload, {
      #   ns <- session$ns
      #   js$collapse(ns("upload_box")) })
      
      # observe({ glimpse(counts()) })
      # observe({ glimpse(metadata()) })
      # observe({ base::colnames(metadata()) })
      # print(metadata_fixed_cols)
      # print(colnames(metadata()))
      
      ##### ** edit metadata #####
      # render metadata table
      output$metadata_table <- renderRHandsontable({
        req(metadata())
        rhandsontable(metadata(), useTypes = TRUE) %>%
          hot_table(highlightCol = TRUE, highlightRow = TRUE, stretchH = "all") %>%
          hot_cols(columnSorting = TRUE) %>%
          hot_col(1:2, readOnly = TRUE)
      })
      
      metadata_final <- reactive({
        if(!is.null(input$metadata_table)) {
          hot_to_r(input$metadata_table)
        } else metadata()
      })
      
      # observe({ glimpse(metadata_final()) })
      
      #### ** edit coldata #####
      # render choose factor UI
      output$choose_factors_ui <- renderUI({
        req(metadata())
        ns <- session$ns
        checkboxGroupInput(ns("choose_factors"), label = NULL, inline = TRUE,
                           choiceNames = colnames(metadata())[-metadata_fixed_cols],
                           choiceValues = colnames(metadata())[-metadata_fixed_cols],
                           selected = colnames(metadata())[-metadata_fixed_cols])
      })
      
      # colData
      coldata <- reactive({
        data <- metadata_final()
        data$condition <- apply(data[, input$choose_factors, drop = FALSE], 1, paste0, collapse = ".")
        data <- filter(data, include == "yes") %>%
          .[, c("sample_name", "batch", "condition"), drop = FALSE] %>%
          column_to_rownames("sample_name")
        return(data)
      })

      # display colData
      output$coldata_table <- renderRHandsontable({
        req(metadata_final())
        rhandsontable(coldata(), readOnly = TRUE, stretchH = "all", rowHeaderWidth = 200) %>%
          hot_cols(columnSorting = TRUE)
      })
      
      # # enable button
      # observe({
      #   req(counts(), metadata_final())
      #   enable("deseq2_apply") })
      
      # observe({ print(glimpse(counts())) })
      # observe({ print(glimpse(metadata_final())) })
      
      # #### ** dds #####
      # dds <- eventReactive(input$deseq2_apply, {
      #   req(counts(), coldata())
      #   counts <- counts()
      #   coldata <- coldata()
      # 
      #   if(!all(rownames(coldata) %in% colnames(counts))) {
      #     showNotification("Not all selected metadata sample_name contained in counts data, check it!", duration = 5, type = "error")
      #     return(NULL)
      #   } else {
      # 
      #     # filter sample according to coldata
      #     counts <- counts[, c(colnames(counts)[1], rownames(coldata))] %>%
      #       column_to_rownames(var=colnames(counts)[1])
      # 
      #     # progress bar
      #     progress <- Progress$new(session, min=1, max=100)
      #     on.exit(progress$close())
      #     progress$set(message = "Construct DESeq2 object", detail = "deseq2 dataset", value = 1)
      # 
      #     # DESeq2
      #     has_batch <- if_else(length(unique(coldata$batch)) > 1, TRUE, FALSE)
      #     # if has batch
      #     if(has_batch) {
      #       deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ batch + condition)
      #     } else {
      #       deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
      #     }
      # 
      #     progress$set(detail = "This takes a while ...", value = 70)
      #     
      #     dds <- DESeq(deseq_dataset, test = "Wald", quiet = TRUE)
      #     
      #     progress$set(detail = "Done", value = 100)
      #     
      #     dds
      #   }
      # })

      # return( reactive({ dds() }) )
      return(list(counts = counts, coldata = coldata))
    }
  )
}
