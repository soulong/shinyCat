

upload <- function(input, output, session, 
                   xlsx_default_sheet = 2 # only work for xlsx file
                   ) {
  
 
  
  return(df)
}



upload_UI <- function(id, 
                      label = NULL,
                      width = 6, 
                      accept_type = c(".csv", ".xlsx")) {
  ns <- NS(id)
  box(
      width = width,
      title = label,
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      fileInput(ns("file"), label = paste0("Choose file (", paste0(accept_type, collapse = "/"), ")"), accept = accept_type),
      fluidRow(
        column(9, uiOutput(ns("select_sheet"))),
        column(3, hidden(actionButton(ns("submit"), label = "Submit")))
      )
  )
}

