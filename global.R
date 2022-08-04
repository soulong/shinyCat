
## load packages
libs <- c("shiny", "shinydashboard", 
          "readxl", "rhandsontable", "rlist", 
          "DT", 
          "shinyjs", 
          "plotly",
          "DESeq2", 
          "corrplot",
          "pheatmap",
          "ggforce",
          "ggrepel",
          "furrr",
          "BiocParallel",
          "clusterProfiler",
          "msigdbr",
          "tidyverse")

invisible(suppressPackageStartupMessages(
  lapply(libs, library, character.only = T)))



## options
options(shiny.maxRequestSize=50*1024^2) # 50M limit size

BiocParallel::register(SnowParam(4)) # multi core
future::plan(multisession, workers = 4)


file_accept_type <- c(".csv", ".xlsx") # expr file extension

counts_data_anno_cols <- 1:4 # 
metadata_fixed_cols <- 1:4 #

deseq2_counts_thred <- 10 # filter rows before deseq2
deseq2_counts_group_number_fraction <- 0.75 # filter rows before deseq2 (group-wise)

digits <- 3 # digits length



# jscode_collapse <- "
# shinyjs.collapse = function(boxid) {
# $('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
# }
# "

