



server <- function(input, output, session) {
  ##### ** upload data #####
  
  # count data
  upload_counts_infile <- reactive({ 
    validate(need(input$counts_data, message = FALSE)) 
    return(input$counts_data) })
  upload_counts_data_filetype <- reactive({
    if_else(last(str_split(upload_counts_infile()$name, "[.]", simplify = TRUE)) == "xlsx", "xlsx", "csv") })
  # select count data sheet UI
  observeEvent(upload_counts_data_filetype(), {
    if(upload_counts_data_filetype() == "xlsx") {
      counts_data_sheets <- excel_sheets(upload_counts_infile()$datapath)
      output$counts_data_select_sheet <- renderUI({
        radioButtons("counts data sheet", label = NULL, choices = counts_data_sheets, selected = last(counts_data_sheets), inline = TRUE) })
    } else {
      output$counts_data_select_sheet <- NULL }
  })
  
  # meta data
  upload_meta_infile <- reactive({ 
    validate(need(input$meta_data, message = FALSE)) 
    return(input$meta_data) })
  upload_meta_data_filetype <- reactive({
    if_else(last(str_split(upload_meta_infile()$name, "[.]", simplify = TRUE)) == "xlsx", "xlsx", "csv") })
  # select meta data sheet UI
  observeEvent(upload_meta_data_filetype(), {
    if(upload_meta_data_filetype() == "xlsx") {
      meta_data_sheets <- excel_sheets(upload_meta_infile()$datapath)
      output$meta_data_select_sheet <- renderUI({
        radioButtons("sample info sheet", label = NULL, choices = meta_data_sheets, selected = first(meta_data_sheets), inline = TRUE) })
    } else {
      output$meta_data_select_sheet <- NULL }
  })
  
  # enable submit button
  observe({ if(!is.null(upload_counts_infile()$datapath)) enable("upload_data_submit") })
  
  # read data
  upload_counts <- eventReactive(input$upload_data_submit, {
    if(is.null(upload_counts_infile()$datapath)) return(NULL)
    showNotification("Read counts data ...", duration = 2, type = "message")
    if(upload_counts_data_filetype() == "xlsx") {
      read_xlsx(upload_counts_infile()$datapath, sheet = input$count_data_select_sheet)
    } else {
      read_csv(upload_counts_infile()$datapath) }
  })
  
  # read meta data
  upload_meta_data <- eventReactive(input$upload_data_submit, {
    if(is.null(upload_meta_infile()$datapath)) return(NULL)
    showNotification("Read meta data ...", duration = 1, type = "message")
    if(upload_meta_data_filetype() == "xlsx") {
      read_xlsx(upload_meta_infile()$datapath, sheet = input$meta_data_select_sheet, col_types = "text") %>%
        mutate(include = as.logical(include))
    } else {
      read_csv(upload_meta_infile()$datapath, col_types = "text")  %>%
        mutate(include = as.logical(include)) }
  })
  
  # collapse upload box
  #  observeEvent(input$upload_data_submit, { js$collapse("upload_box") })
  
  
  
  ##### ** edit meta data #####
  meta_data <- reactive({
    if(!is.null(input$meta_data_table)) {
      hot_to_r(input$meta_data_table)
    } else {
      upload_meta_data() } 
  })
  
  # render table
  output$meta_data_table <- renderRHandsontable({
    rhandsontable(meta_data(), useTypes = TRUE) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, stretchH = "all") %>%
      hot_cols(columnSorting = TRUE) %>%
      hot_col(1:2, readOnly = TRUE) %>%
      hot_context_menu(
        customOpts = list(
          csv = list(name = "Download to CSV",
                     callback = htmlwidgets::JS(
                       "function (key, options) {
                         var csv = csvString(this, sep=',', dec='.');

                         var link = document.createElement('a');
                         link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                         link.setAttribute('download', 'data.csv');

                         document.body.appendChild(link);
                         link.click();
                         document.body.removeChild(link);
                       }"))))
  })
  
  # render choose factor UI
  output$choose_factors_ui <- renderUI({
    checkboxGroupInput("choose_factors", label = NULL, inline = TRUE, 
                       choiceNames = colnames(meta_data())[-meta_data_fixed_cols], 
                       choiceValues = colnames(meta_data())[-meta_data_fixed_cols],
                       selected = colnames(meta_data())[-meta_data_fixed_cols])
  })
  
  # colData
  col_data <- reactive({
    data <- meta_data()
    data$condition <- apply(data[, input$choose_factors], 1, paste0, collapse = ".")
    data <- bind_cols(data[, c(meta_data_fixed_cols)], data[, "condition", drop = FALSE]) %>%
      filter(include == TRUE)
    return(data)
  })
  
  # display colData
  output$col_data_table <- renderRHandsontable({ 
    rhandsontable(col_data(), readOnly = TRUE) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, stretchH = "all") %>%
      hot_cols(columnSorting = TRUE)
  })
  
  # remove samples not include (inculde=no) and convert to matrix
  counts <- eventReactive(input$upload_data_apply, {
    if(!all(col_data()$sample_name %in% colnames(upload_counts()))) {
      showNotification("Not all selected metadata sample_name contained in counts data, check it!", duration = 5, type = "error")
      return(NULL)
    } else {
      upload_counts[, c(colnames(upload_counts)[1], col_data()$sample_name)] %>% 
        column_to_rownames(var="ensembl_gene_id") }
  })
  
  # if has batch
  has_batch <- reactive({ if_else(length(unique(col_data()$batch)) > 1, TRUE, FALSE) })
  
  
  
  
  ##### ** sample qc #####
  
  
}
