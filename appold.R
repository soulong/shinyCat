
library(shiny)
library(shinydashboard)
library(readxl)
library(rlist)
library(rhandsontable)
library(DT)
library(shinyjs)
library(DESeq2)
library(plotly)
library(BiocParallel)
library(tidyverse)
library(clusterProfiler)



##### app options #####
options(shiny.maxRequestSize=50*1024^2) # 50M limit size

register(MulticoreParam(2))

file_accept_type <- c(".csv", ".xlsx")

counts_data_anno_cols <- 1:4
meta_data_fixed_cols <- 1:4

deseq2_counts_thred <- 5
deseq2_counts_group_number_fraction <- 0.75

digits <- 3 # digits length

jscode_collapse <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

# functions
plot_valcano <- function(res, 
                         plot_title="Vocano of DEGs",
                         size=1.5, alpha=0.6,
                         baseMean_thred=10,
                         fc_thred=2, 
                         padj_thred=0.05, 
                         xlims=NULL, ylims=NULL, # or xlims=c(-10, 10), ylims=c(0, 20)
                         legend=FALSE,
                         ns_resampling=500,
                         color=c("blue2", "gray50", "red2")
) {
  df <- filter(res, baseMean > baseMean_thred, !is.na(padj)) %>%
    mutate(`-log10(padj)`=-log10(padj))
  # take sig genes
  df.sig <- filter(df, padj < padj_thred, abs(log2FoldChange) > log2(fc_thred)) %>%
    mutate(direction=if_else(log2FoldChange > 0, "up", "down"))
  # sample un-sig genes,  fo reducing rendering points
  df.unsig <- filter(df, padj >= padj_thred) %>%
    sample_n(size=ns_resampling) %>%
    mutate(direction="ns")
  # merge sig and un-sig
  df.merge <- bind_rows(df.sig, df.unsig)
  if(!is.null(xlims))
    df.merge <- mutate(df.merge, log2FoldChange=if_else(log2FoldChange > xlims[2], xlims[2], 
                                                        if_else(log2FoldChange < xlims[1], xlims[1], log2FoldChange)))
  if(!is.null(ylims))
    df.merge <- mutate(df.merge, `-log10(padj)`=if_else(`-log10(padj)` > ylims[2], ylims[2], 
                                                        if_else(`-log10(padj)` < ylims[1], ylims[1], `-log10(padj)`)))
  # plot
  plot <- ggplot(df.merge, aes(log2FoldChange, `-log10(padj)`, fill=direction,
                               group = 1, text = paste("baseMean: ", baseMean, "</br>external_gene_name: ", external_gene_name))) +
    geom_point(size=size, alpha=alpha, shape = 21, color="transparent") + # shape = 21, color="transparent" will make no stoke
    scale_fill_manual(values=c(down=color[1], ns=color[2], up=color[3])) +
    geom_vline(xintercept=c(-log2(fc_thred), log2(fc_thred)), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thred), linetype="dashed") +
    labs(x="log2 (FoldChange)", y="-log10 (p.adj)",
         title=plot_title, 
         subtitle=str_glue("padj_threshold: ", {{padj_thred}}, 
                           ", FoldChange_threshold: ", {{fc_thred}}, 
                           ", baseMean_threshold: ", {{baseMean_thred}})) +
    cowplot::theme_cowplot()
  if(!is.null(xlims)) plot <- plot + xlim(xlims)
  if(!is.null(ylims)) plot <- plot + ylim(ylims)
  if(legend) plot <- plot + theme(legend.position = "none")
  return(plot)
}

plot_scatter <- function(x) {x}



##### app server #####
server <- function(input, output, session) {

  ##### ** upload data #####
  # count data
  upload_counts_infile <- reactive({
    validate(need(input$counts_data, message = FALSE))
    return(input$counts_data)
  })

  upload_counts_data_filetype <- reactive({
    if_else(last(str_split(upload_counts_infile()$name, "[.]", simplify = TRUE)) == "xlsx", "xlsx", "csv") })
  # select count data sheet UI
  observeEvent(upload_counts_data_filetype(), {
    if(upload_counts_data_filetype() == "xlsx") {
      counts_data_sheets <- excel_sheets(upload_counts_infile()$datapath)
      output$counts_data_select_sheet_ui <- renderUI({
        radioButtons("count_data_select_sheet", label = NULL,
                     choiceNames = counts_data_sheets, choiceValues = counts_data_sheets,
                     selected = dplyr::last(counts_data_sheets), inline = TRUE) })
    } else {
      output$counts_data_select_sheet_ui <- NULL }
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
      output$meta_data_select_sheet_ui <- renderUI({
        radioButtons("meta_data_select_sheet", label = NULL,
                     choiceNames = meta_data_sheets, choiceValues = meta_data_sheets,
                     selected = dplyr::first(meta_data_sheets), inline = TRUE) })
    } else {
      output$meta_data_select_sheet_ui <- NULL }
  })

  # enable submit button if counts was uploaded
  # observe({ if(!is.null(upload_counts_infile()$datapath)) enable("upload_data_submit") })
  observe({
    if(!is.null(upload_counts_infile()$datapath))
      req(upload_counts_infile())
    enable("upload_data_submit")
  })

  # read counts
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

  # view counts
  observeEvent(input$upload_data_submit, { str(upload_counts()) })

  # collapse upload box
  observeEvent(input$upload_data_submit, { js$collapse("upload_box") })


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
    data$condition <- apply(data[, input$choose_factors, drop = FALSE], 1, paste0, collapse = ".")
    data <- filter(data, include == TRUE) %>%
      .[, c("sample_name", "batch", "condition"), drop = FALSE] %>%
      column_to_rownames("sample_name")
    return(data)
  })

  # display colData
  output$col_data_table <- renderRHandsontable({
    rhandsontable(col_data(), readOnly = TRUE, stretchH = "all", rowHeaderWidth = 200) %>%
      hot_cols(columnSorting = TRUE) %>%
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

  # enable run_deseq2_apply buttom
  observe({
    req(upload_counts(), col_data())
    enable("run_deseq2_apply")
  })



  ###### ** run DESeq2 #####
  dds <- eventReactive(input$run_deseq2_apply, {

    col_data <- col_data()
    counts <- upload_counts()

    if(!all(rownames(col_data) %in% colnames(counts))) {
      showNotification("Not all selected metadata sample_name contained in counts data, check it!", duration = 5, type = "error")
      return(NULL)
    } else {
      # filter sample acccording to col_data
      counts <- counts[, c(colnames(counts)[1], rownames(col_data))] %>%
        column_to_rownames(var=colnames(counts)[1])

      # progress bar
      progress <- Progress$new(session, min=1, max=100)
      on.exit(progress$close())
      progress$set(message = "Start construct DESeq2 object", detail = "deseq2 dataset", value = 1)
      Sys.sleep(1)

      # DESeq2
      has_batch <- if_else(length(unique(col_data$batch)) > 1, TRUE, FALSE)
      # if has batch
      if(has_batch) {
        deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ batch + condition)
      } else {
        deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
      }

      # filter low counts
      progress$set(detail = "Filter low counts", value = 30)
      Sys.sleep(1)
      # filter, at least counts_min_sample samples with a count of 10 or higher
      counts_min_sample <- group_by(col_data, condition) %>%
        summarise(n=n()) %>%
        pull(n) %>%
        min()*deseq2_counts_group_number_fraction %>%
        ceiling()
      rows_keep <- rowSums(counts(deseq_dataset) >= deseq2_counts_thred) >= counts_min_sample
      deseq_dataset <- deseq_dataset[rows_keep, ]
      progress$set(detail = paste0("Remaining genes: ", dim(deseq_dataset)[1]), value = 50)
      Sys.sleep(1)

      # run DESeq2
      progress$set(detail = "Run DESeq2, this takes a while ...", value = 70)
      dds <- DESeq(deseq_dataset, test = "Wald", quiet = TRUE)
      progress$set(detail = "Done", value = 100)
      Sys.sleep(0.5)

      # return
      return(dds)
    }
  })




  ###### ** PCA #####
  vst <- eventReactive(input$pca_apply_button, {
    req(dds())
    print(as.logical(input$vst_blind))
    dds()
    vst <- DESeq2::vst(dds(), blind=as.logical(input$vst_blind))
    return(vst)
  })

  output$pca_plot <- renderPlot({
    req(vst())

    p <- DESeq2::plotPCA(vst(), intgroup="condition", ntop=2000, returnData=F) +
      stat_ellipse(aes(fill=group), geom="polygon", alpha=1/4, level=0.9) +
      theme_classic() +
      theme(legend.position="top")
    if(as.logical(input$pca_show_label)) p <- p + ggrepel::geom_text_repel(aes(label=name))

    print(p)
  })


  # enable buttom
  observe({
    req(dds())
    enable("dge_compare_add_buttom")
    enable("dge_compare_submit_buttom")
  })



  ##### ** input comapre UI #####
  dge_compare_counter <- reactiveVal(0)

  # add compare group UI
  observeEvent(input$dge_compare_add_buttom, {
    # add counter, omit 0, start counter from 1
    if(dplyr::first(dge_compare_counter()) == 0) {
      dge_compare_counter(dge_compare_counter() + 1)
    } else {
      dge_compare_counter(c(dge_compare_counter(), last(dge_compare_counter()) + 1))
    }
    n_now <- last(dge_compare_counter())

    # insert UI
    # col_data_conditions <- reactiveVal(c("a", "b", "c")) # for test
    col_data_conditions <- unique(col_data()$condition)
    insertUI(selector = "#dge_compare_add_submit_buttom", where = "beforeBegin", immediate = TRUE,
             ui =
               fluidRow(
                 div(id = paste0("dge_compare_insertUI_", n_now),
                     column(5,
                            selectInput(paste0("dge_compare_contrast_", n_now), NULL, col_data_conditions, selected = col_data_conditions[2])
                     ),
                     column(5,
                            selectInput(paste0("dge_compare_control_", n_now), NULL, col_data_conditions, selected = col_data_conditions[1])
                     ),
                     column(2,
                            actionButton(paste0("remove_dge_compare_buttom_", n_now), "Remove")
                     )
                 )
               )
    )
  })

  # remove compare group UI
  observe({
    # omit if no comparison
    req(min(dge_compare_counter()) > 0)

    # loop see if any remove buttom pressed
    lapply(dge_compare_counter(), function(i) {
      dge_remove_buttom_value <- input[[paste0("remove_dge_compare_buttom_", i)]]
      # req will go on, only if dge_remove_buttom_value value > 0
      # use isTruthy here, req() or need() will interrupte loop
      if(isTruthy(dge_remove_buttom_value)) {
        # print(paste0("remove_dge_compare_buttom_", i))
        removeUI(selector = paste0("#dge_compare_insertUI_", i), immediate = TRUE)
        # delete related counter value
        dge_compare_counter(dge_compare_counter()[-i])
      }
    })
  })




  ##### ** compare result #####
  dge_compare_result_list <- eventReactive(input$dge_compare_submit_buttom, {
    req(dds())

    # get compare group list
    dge_compare_group_list <- lapply(dge_compare_counter(), function(x) {
      return(c(input[[paste0("dge_compare_contrast_", x)]],  input[[paste0("dge_compare_control_", x)]] ))
      })

    # get compare result
    dds <- dds()
    anno <- upload_counts()[, counts_data_anno_cols]
    dge_compare_result_list <- lapply(dge_compare_group_list, function(x) {
      compare_condition <- c("condition", x[1], x[2])
      # unshrink log2FoldChange
      result_unshrink <- results(dds, contrast = compare_condition, cooksCutoff = FALSE, parallel = TRUE)
      # summary(result_unshrink)
      # shrink log2FoldChange
      result_shrink <- lfcShrink(dds, contrast = compare_condition, res = result_unshrink, type = "ashr", parallel = TRUE, quiet = TRUE)
      # bind annotation cols from uoloaded counts data
      result <- left_join(as_tibble(result_unshrink, rownames="ensembl_gene_id")[, -3],
                          as_tibble(result_shrink, rownames="ensembl_gene_id")[, c(1, 3)],
                          by=c("ensembl_gene_id"="ensembl_gene_id")) %>%
        # transmute(ensembl_gene_id=ensembl_gene_id,
        #           baseMean=round(baseMean), log2FoldChange=round(log2FoldChange, digits),
        #           stat=round(stat, digits), pvalue=signif(stat, digits), padj=signif(padj, digits)) %>%
        dplyr::select(ensembl_gene_id, baseMean, log2FoldChange, stat, pvalue, padj) %>%
        right_join(anno, ., by=c("ensembl_gene_id"="ensembl_gene_id")) %>%
        arrange(pvalue)
      return(result)
    })
    # names compare result by "contrast_vs_control"
    names(dge_compare_result_list) <- lapply(dge_compare_group_list, paste0, collapse = "_vs_") %>% unlist()

    return(dge_compare_result_list)
  })



  ###### ** sig DGEs #####
  dge_compare_result_sig_list <- reactive({
    req(dge_compare_result_list())
    dge_compare_result_sig_list <- lapply(dge_compare_result_list(), function(x) {
      dplyr::filter(x,
                    baseMean > input$dge_baseMean_threshold,
                    abs(log2FoldChange) > log2(input$dge_fc_threshold),
                    padj < input$dge_p.adj_threshold)
    })
    names(dge_compare_result_sig_list) <- names(dge_compare_result_list())
    return(dge_compare_result_sig_list)
  })

  # collapse setup box
  observeEvent(input$dge_compare_submit_buttom, {
    js$collapse("dge_analysis_threshold_setup")
    js$collapse("dge_analysis_add_comparison")
  })



  ##### ** plot DGEs #####
  # DT table, volcano plot, scatter plot
  output$dge_sig_tables_ui <- renderUI({
    validate(need(dge_compare_result_sig_list(), message = FALSE))

    sig_list <- dge_compare_result_sig_list()
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
    do.call(what = shinydashboard::tabBox, args = c(id = "dge_sig_table_tabbox", width = 12, tabpanel_list))
  })

  # observe({
  #   try(print(input$dge_sig_table_tabbox))
  #   try(glimpse(dge_compare_result_list()[[input$dge_sig_table_tabbox]]))
  # })
  #

  # dispaly valcano
  output$dge_sig_valcano <- renderPlotly({
    req(dge_compare_result_list())
    # plotly volcano
    ggplotly(
      plot_valcano(dge_compare_result_list()[[input$dge_sig_table_tabbox]],
                   plot_title = input$dge_sig_table_tabbox,
                   baseMean_thred = input$dge_baseMean_threshold,
                   fc_thred = input$dge_fc_threshold,
                   padj_thred = input$dge_p.adj_threshold),
      tooltip = c("text", "log2FoldChange", "-log10(padj)"))
  })


  ##### ** ORA #####
  # collapge box
  observeEvent(input$ora_apply_enrich, {
    js$collapse("ora_submit_box")
  })

  # ora enrich
  dge_compare_result_enrich_ora_list <- eventReactive(input$ora_apply_enrich, {
    req(dge_compare_result_sig_list())
    orgdb <- if_else(input$ora_species=="hs", "org.Hs.eg.db", "org.Mm.eg.db")

    # progress bar
    progress <- Progress$new(session, min=0, max=2*length(dge_compare_result_sig_list()))
    on.exit(progress$close())
    progress$set(message = "Gene Ontology ORA enrichment ... ", value = 0)

    enrich_list <- lapply(seq_along(dge_compare_result_sig_list()), function(i) {
      up <- filter(dge_compare_result_sig_list()[[i]], log2FoldChange > 0) %>% pull(entrezgene_id)
      down <- filter(dge_compare_result_sig_list()[[i]], log2FoldChange < 0) %>% pull(entrezgene_id)

      progress$set(detail = paste0(names(dge_compare_result_sig_list())[i], " : Up-regulated"), value = (2 * i - 1))
      up_enrich <- enrichGO(up, orgdb, ont = "ALL", qvalueCutoff=0.9,
                            pvalueCutoff = input$ora_pval_cutoff, minGSSize= input$ora_min_setSize, maxGSSize = input$ora_max_setSize, readable = TRUE)
      up_enrich@result$direction <- "Up-regulated"

      progress$set(detail = paste0(names(dge_compare_result_sig_list())[i], " : Down-regulated"), value = (2 * i))
      down_enrich <- enrichGO(down, orgdb, ont = "ALL", qvalueCutoff=0.9,
                            pvalueCutoff = input$ora_pval_cutoff, minGSSize= input$ora_min_setSize, maxGSSize = input$ora_max_setSize, readable = TRUE)
      down_enrich@result$direction <- "Down-regulated"

      result <- bind_rows(up_enrich@result, down_enrich@result)
      return(result)
    })
    names(enrich_list) <- names(dge_compare_result_sig_list())
    return(enrich_list)
  })

  # observe({
  #   print(names(dge_compare_result_enrich_ora_list())[[1]])
  # })

  # display ora enrich
  output$enrich_ora_table_ui <- renderUI({
    validate(need(dge_compare_result_enrich_ora_list(), message = FALSE))

    ora_list <- dge_compare_result_enrich_ora_list()
    tabpanel_list <- lapply(seq_along(ora_list), function(i) {
      tabPanel(title = names(ora_list)[i],
               renderDT(ora_list[[i]], server = FALSE, # server = FALSE, it will download all data, not only visable data
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
    do.call(what = shinydashboard::tabBox, args = c(id = "enrich_ora_table_tabbox", width = 12, tabpanel_list))
  })

  # display ora enrich barplot
  output$enrich_ora_barplot <- renderPlot({
    req(dge_compare_result_enrich_ora_list(), input$enrich_ora_table_tabbox)

    data <- dge_compare_result_enrich_ora_list()[[input$enrich_ora_table_tabbox]] %>%
      group_by(direction, ONTOLOGY) %>%
      arrange(pvalue) %>%
      dplyr::slice(1:10) %>%
      mutate(`-log10(p.adjust)`=-log10(p.adjust)) %>%
      ungroup()
    # dotplot
    ggpubr::ggbarplot(data, x="Description", y="-log10(p.adjust)", fill="ONTOLOGY", color="white",
              orientation="horizontal", sort.val="asc", facet.by="direction", scales="free_y") +
      scale_x_discrete(labels=function(x) str_wrap(x, width=40))
  })



  
  ##### ** GSEA #####
  
  
  


 
}





##### app ui #####
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
      menuItem("DGE analysis", tabName = "dge_analysis")
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
   # extendShinyjs(text = jscode_collapse),
    
    tabItems(
      ##### **** tab upload #####
      tabItem(tabName = "upload",
              fluidRow(
                box(id = "upload_box", title = "Upload data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    fileInput("counts_data",  paste0("Gene expression raw counts (", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
                    uiOutput("counts_data_select_sheet_ui"),
                    fileInput("meta_data", paste0("Sample info metadata (optional, ", paste0(file_accept_type, collapse = "/"), ")"), accept = file_accept_type),
                    uiOutput("meta_data_select_sheet_ui"),
                    disabled(actionButton("upload_data_submit", "Submit"))
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
                             div(align = "center", disabled(actionButton("run_deseq2_apply", "Run"))))
                    ),
                    rHandsontableOutput("col_data_table")
                )
              )
      ),
      
      
      ##### **** tab sample_qc #####
      tabItem(tabName = "sample_qc",
              box(width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE,
                  column(4,
                         radioButtons("vst_blind", "vst transformation blind", selected = "true", inline = TRUE, choices = c(true="true", false="false"))
                  ),
                  column(4,
                         radioButtons("pca_show_label", "show labels on plot", selected = "true", inline = TRUE, choices = c(true="true", false="false"))
                  ),
                  column(4, align="center",
                         actionButton("pca_apply_button", "Plot PCA")
                  )
              ),
              box(width = 12, title = "sample PCA", 
                  plotOutput("pca_plot")
              )
      ),
      
      
      ##### **** tab dge_analysis #####
      tabItem(tabName = "dge_analysis",
              fluidRow(
                column(7, 
                       box(id = "dge_analysis_add_comparison", title = "Choose comparisons", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                           div(id = "dge_compare_add_submit_buttom", align = "center", 
                               column(6, align = "center", 
                                      disabled(actionButton("dge_compare_add_buttom", "Add new comparison"))),
                               column(6, align = "center", 
                                      disabled(actionButton("dge_compare_submit_buttom", "Get compare result")))
                           )
                       )
                ),
                column(5, 
                       box(id = "dge_analysis_threshold_setup", title = "Significance threshold", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                           column(4, align = "left", 
                                  numericInput("dge_baseMean_threshold", "baseMean", 1, 1, Inf, width = "100%")
                           ),
                           column(4, align = "left", 
                                  numericInput("dge_fc_threshold", "FoldChange", 1.5, 1, Inf, width = "100%")
                           ),
                           column(4, align = "left",
                                  numericInput("dge_p.adj_threshold", "Adjust.pvalue", 0.05, 0, Inf, width = "100%")
                           )
                       )
                )
              ),
              fluidRow(
                uiOutput("dge_sig_tables_ui")
              ),
              fluidRow(   
                box(width = 12, 
                    column(6, 
                           plotlyOutput("dge_sig_valcano")
                    ),
                    column(6,
                           
                    )
                )
              )
      ),

      ##### **** tab ora #####
      tabItem(tabName = "ora",
              box(id = "ora_submit_box", title = "ORA enrichment setup", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                  column(2, align="center", radioButtons("ora_species", "Species", selected ="mm", inline = TRUE, choiceNames = c("hs", "mm"), choiceValues = c("hs", "mm"))
                  ),
                  column(2, numericInput("ora_pval_cutoff", "ora_pval_cutoff", 0.1, min = 0, max = 0.2)
                  ),
                  column(2, numericInput("ora_min_setSize", "ora_min_setSize", 10, min = 5, max = 50)
                  ),
                  column(2, numericInput("ora_max_setSize", "ora_max_setSize", 400, min = 50, max = 1000)
                  ),
                  column(4, align="center", br(), actionButton("ora_apply_enrich", "Enrich now", width = "200px")
                  )
              ),
              fluidRow(
                uiOutput("enrich_ora_table_ui")
              ),
              fluidRow(
                plotOutput("enrich_ora_barplot")
              )
      ),
      
      ##### **** tab gsea #####
      tabItem(tabName = "gsea",
              fluidRow(
                
              )
      )
      
    ) # end of tabItems
  ) # end of dashboardBody
  
) # end of dashboardPage



##### runapp #####
# shinyApp(ui = ui, server = server)
