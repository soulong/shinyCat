
# module_ui <- function(id, xyz) {
#   ns <- NS(id)
#   tagList(
#     fileInput(ns("id"), zyx),
#     ...
#   )
# }
sample_qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 6, title = "DESeq2 option", status = "warning", solidHeader = TRUE, collapsible = TRUE,
          column(4,
                 radioButtons(ns("filter_symbol"), "Remove rows with empty symbol(gene_name)", selected = "true", inline = TRUE, choices = c(true="true", false="false"))
          ),
          column(4,
                 numericInput(ns("min_counts"), "Counts threshold", value = 10, min = 1, step = 1)
          ),
          column(4, align = "center",
                 actionButton(ns("apply"), "Apply"),
          )
      ),
      box(width = 6, title = "Plot option", status = "warning", solidHeader = TRUE, collapsible = TRUE,
          column(6,
                 radioButtons(ns("pca_show_label"), "Labels on PCA", selected = "true", inline = TRUE, choices = c(true="true", false="false"))
          ),
          column(6,
                 radioButtons(ns("calculate_source"), "Calculating use", selected = "vst", inline = TRUE, choices = c(vst="vst", counts_norm="counts_norm"))
          )
      )
    ),
    fluidRow(
      box(width = 6, title = "PCA & Correlation", status = "primary", solidHeader = TRUE,
          plotOutput(ns("pca_plot")),
          br(),
          plotOutput(ns("cor_plot"))
      ),
      box(width = 6, title = "Top 2000 variable genes", status = "primary", solidHeader = TRUE,
          plotOutput(ns("top_var_plot"), height = "800px"),
          br()
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
sample_qc_server <- function(id, counts, coldata) {
  moduleServer(
    id,
    function(input, output, session) {
      
      #### ** dds #####
      dds <- eventReactive(input$apply, {
        req(counts(), coldata())
        counts <- counts()
        coldata <- coldata()
        
        if(!all(rownames(coldata) %in% colnames(counts))) {
          showNotification("Not all selected metadata sample_name contained in counts data, check it!", duration = 5, type = "error")
          return(NULL)
        } else {
          
          # progress bar
          progress <- Progress$new(session, min=1, max=10)
          on.exit(progress$close())

          
          # filter gene_name
          progress$set(message = "Filter empty symbols", value = 2)
          if(as.logical(input$filter_symbol)) {
            counts <- filter(counts, !is.na(gene_name))
          } else {
            counts <- mutate(counts, gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
          }
          
          # filter sample according to coldata
          progress$set(message = "Filter samples acorrding to colData", value = 3)
          counts <- counts[, c("gene_name", rownames(coldata))] %>%
            distinct(gene_name, .keep_all = T) %>%
            column_to_rownames(var = "gene_name")
          
          # filter counts
          progress$set(message = "Filter low counts rows", value = 4)
          counts <- counts[rowMeans(counts) > input$min_counts, ]
          # print(counts[1:5,1:5])
          
          # DESeq2
          progress$set(message = "Construct DESeq2 object",  value = 5)
          if(length(unique(coldata$batch)) > 1) {
            # if has batch
            deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ batch + condition)
          } else {
            deseq_dataset <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
          }
          
          progress$set(detail = "This takes a while ...", value = 6)
          
          dds <- DESeq(deseq_dataset, test = "Wald", quiet = TRUE)
          
          progress$set(detail = "Done", value = 100)
          
          dds
        }
      })
      
      mat <- reactive({
        req(dds())
        d <- dds()
        if(input$calculate_source == "vst") {
          assay(DESeq2::vst(d, blind = FALSE))
        } else {
          DESeq2::counts(d, normalized = TRUE)
        }
      })
      
      # observe({ mat()[1:5,1:5] })
      
      ##### ** PCA ##### 
      output$pca_plot <- renderPlot({
        req(dds())
        pcaData <- plotPCA(vst(dds()), intgroup=c("condition"), returnData=TRUE)
        # glimpse(pcaData)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
          geom_point(size=2) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
          geom_mark_ellipse(aes(fill = condition, color = condition), expand = unit(1, "mm")) +
          # stat_ellipse(aes(fill=condition), geom="polygon", alpha=0.2, level=0.9) +
          coord_fixed() +
          theme_bw(16) +
          theme(legend.position="top")
        # shouw labels
        if(as.logical(input$pca_show_label)) p <- p + geom_text_repel(aes(label = name))
        print(p)
      })
      
      ###### ** correlation #####
      output$cor_plot <- renderPlot({
        req(mat())
        corrplot.mixed(cor(mat()), 
                       lower = "number", upper="ellipse", 
                       order ="AOE", 
                       tl.pos = "lt", tl.col="black",
                       number.cex=0.7)
      })
      
      ###### ** variable genes heatmap #####
      output$top_var_plot <- renderPlot({
        req(mat())
        top_variable_genes <- order(rowVars(mat()), decreasing = TRUE)[1:2000]
        df <- as.data.frame(colData(dds())[, "condition", drop = FALSE])
        pheatmap(mat()[top_variable_genes, ], 
                 scale = "row",
                 show_rownames=FALSE,
                 fontsize_col = 12,
                 angle_col = 45,
                 cluster_cols=FALSE, 
                 annotation_col=df)
        })

      
      # observe({ write_rds(dds(), "dds.rds") })
      return(dds)
    }
  )
}
