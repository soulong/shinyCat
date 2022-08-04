
# module_ui <- function(id, xyz) {
#   ns <- NS(id)
#   tagList(
#     fileInput(ns("id"), zyx),
#     ...
#   )
# }
deseq2_compare_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(7, 
             box(id = ns("deg_compare"), title = "Choose comparisons", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 actionButton(ns("add_compare"), "Add new comparison"),
                 actionButton(ns("apply"), "Get compare result")
             )
      ),
      column(5, 
             box(title = "Significance threshold", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 column(4, align = "left", 
                        numericInput(ns("baseMean_threshold"), "baseMean", 100, 1, Inf, width = "100%")
                 ),
                 column(4, align = "left", 
                        numericInput(ns("fc_threshold"), "FoldChange", 1.5, 1, Inf, width = "100%")
                 ),
                 column(4, align = "left",
                        numericInput(ns("padj_threshold"), "Adjust.pvalue", 0.05, 0, Inf, width = "100%")
                 )
             )
      )
    ),
    fluidRow(
      uiOutput(ns("sig_tables_ui"))
    ),
    fluidRow(   
      box(width = 12, 
          column(6, 
                 plotlyOutput(ns("sig_valcano"))
          ),
          column(6,
                 
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
deseq2_compare_server <- function(id, dds) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ##### ** input comapre UI #####
      counter <- reactiveVal(0)
      
      # add compare group UI
      observeEvent(input$add_compare, {
        req(dds())
        ns = session$ns
        
        # add counter, omit 0, start counter from 1
        if(first(counter()) == 0) {
          counter(counter() + 1)
        } else counter(c(counter(), last(counter()) + 1))
        n_now <- last(counter())
        # print(counter())
        
        # insert UI
        # conditions <- reactiveVal(c("a", "b", "c")) # for test
        conditions <- unique(colData(dds())$condition)
        insertUI(selector = paste0("#", ns("add_compare")), where = "beforeBegin", immediate = TRUE,
                 ui =
                   fluidRow(
                     div(id = ns(paste0("compare_insertUI_", n_now)),
                         column(5,
                                selectInput(ns(paste0("compare_contrast_", n_now)), NULL, conditions, selected = conditions[2])
                         ),
                         column(5,
                                selectInput(ns(paste0("compare_control_", n_now)), NULL, conditions, selected = conditions[1])
                         ),
                         column(2,
                                actionButton(ns(paste0("remove_compare_", n_now)), "Remove")
                         )
                     )
                   )
        )
      })
      

      # remove compare group UI
      observe({
        # omit if no comparison
        req(min(counter()) > 0)
        ns = session$ns

        # loop see if any remove buttom pressed
        lapply(counter(), function(i) {
          remove_compare_buttom <- input[[paste0("remove_compare_", i)]]
          # req will go on, only if deg_remove_buttom_value value > 0
          # use isTruthy here, req() or need() will interrupt loop
          if(isTruthy(remove_compare_buttom)) {
            # print(paste0("remove_deg_compare_buttom_", i))
            removeUI(selector = paste0("#", ns(paste0("compare_insertUI_", i))), immediate = TRUE)
            # delete related counter value
            counter(counter()[-i])
          }
          # print(counter())
        })
      })
      
      
      ##### ** compare result #####
      results <- eventReactive(input$apply, {
        req(min(counter()) > 0)
        d <- dds()
        
        # progress bar
        progress <- Progress$new(session, min=1, max=10)
        on.exit(progress$close())
        progress$set(message = "Wald test", value = 1)
        
        # get compare group list
        compare_list <- map(counter(), function(i) {
          c(input[[paste0("compare_contrast_", i)]], input[[paste0("compare_control_", i)]] ) })
        
        # get compare result
        # anno <- upload_counts()[, counts_data_anno_cols]
        results <- map(compare_list, function(i) {
          compare_condition <- c("condition", i[1], i[2])
          progress$set(detail = compare_condition, value = 5)
          
          result_unshrink <- DESeq2::results(d, contrast = compare_condition)
          result_shrink <- lfcShrink(d, contrast = compare_condition, res = result_unshrink, type = "ashr", quiet = TRUE)

          # tidy
          res <- as_tibble(result_unshrink, rownames = "gene_name")[, -3] %>%
            left_join(as_tibble(result_shrink, rownames = "gene_name")[, c(1, 3)], by = c("gene_name"="gene_name")) %>%
            relocate(log2FoldChange, .before = pvalue) %>%
            arrange(pvalue) %>%
            mutate(baseMean = round(baseMean, 0), 
                   log2FoldChange = round(log2FoldChange, 5), 
                   lfcSE = round(lfcSE, 5),
                   stat = round(stat, 5))
          # glimpse(res)
          return(res)
        })
        
        progress$set(detail = "Done", value = 10)

        names(results) <- map(compare_list, paste0, collapse = "_vs_") %>% 
          unlist()
        
        return(results)
      })
      
      # observe({
      #   glimpse(results) })
      
      
      ###### ** filter significant #####
      results_sig <- reactive({
        req(results())

        results_sig <- map(
          results(), 
          ~ dplyr::filter(.x, 
                          baseMean > input$baseMean_threshold,
                          abs(log2FoldChange) > log2(input$fc_threshold),
                          padj < input$padj_threshold
          ))
        # names(results_sig) <- names(results())
        return(results_sig)
      })

      # observe({
      #   glimpse(results_sig) })

      # # collapse setup box
      # observeEvent(input$apply, {
      #   js$collapse("deg_analysis_threshold_setup")
      #   js$collapse("deg_analysis_add_comparison")
      # })
      
      
      ##### ** results #####
      # DT table, volcano plot, scatter plot
      output$sig_tables_ui <- renderUI({
        req(results_sig())
        ns = session$ns

        sig_list <- results_sig()
        tabpanel_list <- lapply(seq_along(sig_list), function(i) {
          tabPanel(title = names(sig_list)[i],
                   renderDT(sig_list[[i]], server = FALSE, # server = FALSE, it will download all data, not only visable data
                            filter = "top", selection = "single", rownames = FALSE,
                            style = "bootstrap", class = "cell-border stripe",
                            extensions = c("Buttons", "Scroller"),
                            options = list(dom = "Brtip", # Bfrtip will display searching bar
                                           columnDefs = list(list(className = "dt-center", targets = "_all")),
                                           buttons = c(I("colvis"), "copy", "excel"), # Buttons extension
                                           deferRender = TRUE, scrollY = "150px", scroller = TRUE, # Scroller extension
                                           scrollX = TRUE, searching = TRUE)
                   ))
        })
        do.call(what = shinydashboard::tabBox, args = c(id = ns("sig_table_tabbox"), width = 12, tabpanel_list))
      })

      # observe({
      #   try(print(input$deg_sig_table_tabbox))
      #   try(glimpse(deg_compare_result_list()[[input$deg_sig_table_tabbox]]))
      # })
      #

      ##### ** valcano #####
      output$sig_valcano <- renderPlotly({
        req(results_sig(), input$sig_table_tabbox)

        df <- results_sig()[[input$sig_table_tabbox]]
        # glimpse(df)
        ggplotly(
          plot_valcano(df,
                       plot_title = input$sig_table_tabbox,
                       baseMean_thred = input$baseMean_threshold,
                       fc_thred = input$fc_threshold,
                       padj_thred = input$padj_threshold),
          tooltip = c("text", "log2FoldChange", "-log10(padj)"))
      })

      observe({ write_rds(list(results = results(), results_sig = results_sig()), "results.rds") })
      return(list(results = results, results_sig = results_sig))
    }
  )
}
