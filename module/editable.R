

editable <- function(input, output, session, 
                     df, 
                     filter = "none", # top
                     editable = FALSE, # list(target = "cell", disable = list(columns = c(1, 2)))
                     pageLength = 20
                     ) {
  
  # render factor choose UI
  output$choose_factor <- renderUI({
    ns <- session$ns
    checkboxGroupInput(ns("factor"), label = "Choose factor", inline = TRUE,
                       choiceNames = colnames(df())[-1:-4], 
                       choiceValues = colnames(df())[-1:-4])
  })
  
  # colData
  re <- reactive({
    data <- df()
    data$condition <- apply(data[, input$factor], 1, paste, collapse=".")

    # batch
    has_batch <- if_else(unique(data$batch) > 1, TRUE, FALSE)
    # coldata
    if(has_batch) {
      col_data <- dplyr::filter(data, isTRUE(include)) %>%
        dplyr::select(sample_name, condition, batch) %>%
        column_to_rownames("sample_name")
    } else {
      col_data <- dplyr::filter(data) %>%
        dplyr::select(sample_name, condition) %>%
        column_to_rownames("sample_name")
    }
    
    return(list(col_data=col_data, has_batch=has_batch))
  })
  
  
  # render table
  output$table <- renderDT(
    df(), server = TRUE, filter = filter, editable = editable, # scrollX = TRUE,
    extensions = "Buttons",
    options = list(pageLength = pageLength, autoWidth = TRUE, searchHighlight = TRUE,
                   dom = "Bfrtip", buttons = c("copy", "excel")
                   )
    )
  
  # render condtion
  output$conditon_table <- renderDT(
    re()$col_data, rownames = TRUE, filter = "none", 
    options = list(pageLength = pageLength)
  )
  
  # # update new table
  # observeEvent(input$table_cell_edit, {
  #   table_new <- df()
  #   info = input$table_cell_edit
  #   str(info)
  #   i = info$row
  #   j = info$col
  #   glimpse(df())
  #   glimpse(table_new)
  #   table_new[i, j] <- DT::coerceValue(info$value, table_new[i, j])
  #   glimpse(table_new)
  # })

  return(re)
}



editable_UI <- function(id,
                        width = 6, 
                        button_label = "Apply"
                        ) {
  ns <- NS(id)
  tagList(
    box(
      width = width, 
      title = "Metadata",
      status = "warning",
      solidHeader = TRUE,
      collapsible = TRUE,
      # fluidRow(align = "center", actionButton(ns("submit"), label = button_label), width = "400px"),
      div(DTOutput(ns("table")), style= "overflow-x: auto")
    ),
    box(
      width = width, 
      title = "Condition",
      status = "warning",
      solidHeader = TRUE,
      collapsible = TRUE,
      uiOutput(ns("choose_factor")),
      div(DTOutput(ns("conditon_table")), style= "overflow-x: auto")
    )
  )
}
