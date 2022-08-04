
## functions
file_extension <- function(x) str_split(x, "[.]", simplify = TRUE) %>% last()


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
    slice_sample(n = ns_resampling) %>%
    mutate(direction = "ns")
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
                               group = 1, text = paste("baseMean: ", baseMean, "</br>gene_name: ", gene_name))) +
    geom_point(size=size, alpha=alpha, shape = 21, color="transparent", show.legend = legend) + # shape = 21, color="transparent" will make no stoke
    scale_fill_manual(values=c(down=color[1], ns=color[2], up=color[3])) +
    geom_vline(xintercept=c(-log2(fc_thred), log2(fc_thred)), linetype="dashed", color = "grey50") +
    geom_hline(yintercept=-log10(padj_thred), linetype="dashed", color = "grey50") +
    labs(x="log2(FoldChange)", y="-log10(padj)",
         title=plot_title,
         subtitle=str_glue("padj_threshold: ", {{padj_thred}},
                           ", FoldChange_threshold: ", {{fc_thred}},
                           ", baseMean_threshold: ", {{baseMean_thred}})) +
    theme_bw(14)
  if(!is.null(xlims)) plot <- plot + xlim(xlims)
  if(!is.null(ylims)) plot <- plot + ylim(ylims)
  
  return(plot)
}
