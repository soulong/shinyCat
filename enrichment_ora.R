
enrichment_ora_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(title = "ORA enrichment setup", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
        column(2, align="center", radioButtons(ns("species"), "Species", selected ="hs", inline = TRUE, choiceNames = c("hs", "mm"), choiceValues = c("hs", "mm"))
        ),
        column(2, align="center", radioButtons(ns("database"), "database", selected ="HALLMARK", inline = TRUE, choiceNames = c("HALLMARK", "GOBP", "KEGG"), choiceValues = c("HALLMARK", "GOBP", "KEGG"))
        ),
        column(2, numericInput(ns("pval_cutoff"), "ora_pval_cutoff", 0.1, min = 0, max = 0.2)
        ),
        column(2, numericInput(ns("min_setSize"), "ora_min_setSize", 10, min = 5, max = 50)
        ),
        column(2, numericInput(ns("max_setSize"), "ora_max_setSize", 400, min = 50, max = 1000)
        ),
        column(2, align="center", br(), actionButton(ns("apply"), "Enrich", width = "200px")
        )
    ),
    fluidRow(
      uiOutput(ns("enrich_table_ui"))
    ),
    fluidRow(
      plotOutput(ns("enrich_barplot"))
    )
  )
}




enrichment_ora_server <- function(id, results_sig) {
  moduleServer(
    id,
    function(input, output, session) {
      
     
      # db
      genesets <- eventReactive(input$apply, {
        req(results_sig())
        glimpse(results_sig())
        
        species <- if_else(input$species=="hs", "human", "mouse")
        msigdb <- msigdbr::msigdbr(species = species)
        genesets <- case_when(
          input$database == "HALLMARK" ~ filter(msigdb, gs_cat == "H"),
          input$database == "KEGG" ~ filter(msigdb, gs_subcat == "CP:KEGG"),
          input$database == "GOBP" ~ filter(msigdb, gs_subcat == "GO:BP")
        )
        
        # genesets <- list(
        #   HALLMARK = filter(msigdb, gs_cat == "H")
        #   # KEGG = filter(msigdb, gs_subcat == "CP:KEGG"),
        #   # GO = filter(msigdb, gs_subcat == "GO:BP")
        #   # BIOCARTA = filter(msigdb, gs_subcat == "CP:BIOCARTA"),
        #   # PID = filter(msigdb, gs_subcat == "CP:PID")
        # )
        genesets <- map(genesets, 
                        ~ dplyr::select(.x, gs_name, gene_symbol) %>%
                          distinct())
        glimpse(genesets)
        write_rds(genesets, "genesets.rds")
        genesets
      })
      
      observe({ glimpse(genesets()) })
      
      # enrichment
      enrich_list <- eventReactive(input$apply, {
        req(results_sig(), genesets())
        
        # progress bar
        progress <- Progress$new(session, min=0, max=length(results_sig()))
        on.exit(progress$close())
        progress$set(message = "ORA enrichment", value = 0)
        
        enrich_list <- map(seq_along(results_sig()), function(i) {
          # up
          print(i)
          progress$set(detail = paste0(names(results_sig())[i], " : upregulated"), value = i)
          up <- filter(results_sig()[[i]], log2FoldChange > 0) %>% pull(gene_name)
          ora_up <- enricher(up, pvalueCutoff = input$pval_cutoff, qvalueCutoff = 0.7, 
                             minGSSize = input$min_setSize, maxGSSize = input$max_setSize, 
                             TERM2GENE = genesets())
          ora_up@result$direction <- "Up"
          # down
          progress$set(detail = paste0(names(results_sig())[i], " : downregulated"))
          down <- filter(results_sig()[[i]], log2FoldChange > 0) %>% pull(gene_name)
          ora_down <- enricher(down, pvalueCutoff = input$pval_cutoff, qvalueCutoff = 0.7, 
                               minGSSize = input$min_setSize, maxGSSize = input$max_setSize, 
                               TERM2GENE = genesets())
          ora_down@result$direction <- "Down"
          
          # merge
          ora <- bind_rows(ora_up@result, ora_down@result) %>%
            as_tibble()
          
          return(ora)
        })
        
        # names(enrich_list) <- names(results_sig())
        return(enrich_list)
      })
      
      observe({ glimpse(enrich_list()) })
      

      # display ora enrich
      output$enrich_table_ui <- renderUI({
        req(enrich_list())
        ns = session$ns
        
        # ora_list <- enrich_list()
        tabpanel_list <- lapply(seq_along(enrich_list()), function(i) {
          tabPanel(title = names(enrich_list())[i],
                   renderDT(enrich_list()[[i]], server = FALSE, # server = FALSE, it will download all data, not only visable data
                            filter = "top", selection = "single", rownames = FALSE,
                            style = "bootstrap", class = "cell-border stripe",
                            extensions = c("Buttons", "Scroller"),
                            options = list(dom = "Brtip", # Bfrtip will display searching bar
                                           columnDefs = list(list(className = "dt-center", targets = "_all")),
                                           buttons = c(I("colvis"), "copy", "excel"), # Buttons extension
                                           deferRender = TRUE, scrollY = "100px", scroller = TRUE, # Scroller extension
                                           scrollX = TRUE, searching = TRUE)
                   ))
        })
        do.call(what = shinydashboard::tabBox, args = c(id = ns("enrich_table_tabbox"), width = 12, tabpanel_list))
      })
      
      # display ora enrich barplot
      output$enrich_barplot <- renderPlot({
        req(enrich_list(), input$enrich_table_tabbox)
        
        data <- results_sig()[[input$enrich_table_tabbox]] %>%
          group_by(direction) %>%
          arrange(pvalue) %>%
          dplyr::slice(1:10) %>%
          mutate(`-log10(p.adjust)`=-log10(p.adjust)) %>%
          ungroup()
        # dotplot
        ggpubr::ggbarplot(data, x="Description", y="-log10(p.adjust)", fill="ONTOLOGY", color="white",
                          orientation="horizontal", sort.val="asc", facet.by="direction", scales="free_y") +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 60))
      })
     })
}
