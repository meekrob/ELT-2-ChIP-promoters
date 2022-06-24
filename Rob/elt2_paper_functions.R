
# myPDFplot(), wrapper to save non-ggplot images to a PDF file

myPDFplot <- function(plot, name, height, width, plotdir = plotdir) {
  pdf(
    paste(plotdir,
          name,
          "_",
          lubridate::today(),
          ".pdf",
          sep = ""),
    height = height,
    width = width
  )
  print(plot)
  dev.off()
}

# vsd.corr.per.stage(), plot correlation matrix for a given stage
vsd.corr.per.stage <- function(x, main){
  vsd <- assay(vsd)[,metadata1 %>% filter(grepl(x, names)) %>% pull(names)]
  sampleDists <- dist(t(vsd))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors, 
           main = main)
}

# res_to_df(), turn DESeq2 result object to dataframe
res_to_df <- function(in_res){
  out_df <- as.data.frame(in_res) 
  if (! 'WBGeneID' %in% colnames(out_df))
    {out_df %<>% rownames_to_column(var = "WBGeneID")}
  out_df
}

shrunk_pairwise_array_df <- function(stage) {
  samples <- paste(stage, c("whole", "GFPplus", "cells", "GFPminus"), sep = "")
  combos <- combn(samples, 2, simplify = FALSE)
  # print(combos)
  all_pairwise_comparisons <- data.frame()
  for (i in 1:length(combos)) {
    tobind <-
      as.data.frame(lfcShrink(
        dds,
        contrast = c("group", combos[[i]][1], combos[[i]][2]),
        type = "ashr",
        quiet = TRUE
      )) %>%
      rownames_to_column(var = "WBGeneID") %>%
      mutate(comparison = paste(combos[[i]][1], combos[[i]][2], sep = "_vs_"))  %>% 
      mutate(label = str_remove_all(comparison, "embryo|L1|L3"))
    all_pairwise_comparisons <- bind_rows(all_pairwise_comparisons, tobind)
  }
  all_pairwise_comparisons
}

MA_plot_array <- function(in.df, title, sig){
  ggplot(in.df %>% mutate(padj = replace_na(padj, 1)), aes(x = log10(baseMean), y = log2FoldChange, color = padj < sig)) +
    geom_point(shape = 16, alpha = 0.1, stroke = 0, size = 1) +
    ylim(c(-10,10))+
    facet_wrap(~label) +
    scale_color_manual(values = c("black", "red"), name = "q.value < 0.1") +
    theme_classic() +
    ggtitle(title)
}

alt_hyp_res_df <- function(stage, thresh, sig){
  samples <- paste(stage, c("whole", "GFPplus", "cells", "GFPminus"), sep = "")
  combos <- combn(samples, 2, simplify = FALSE)
  hyps = c("greater", "less", "lessAbs")
  df <- data.frame()
  for(i in 1:length(combos)){
    for(hyp in hyps){
      thresh_res <- results(dds, contrast = c("group", combos[[i]][1],combos[[i]][2]), lfcThreshold=thresh, altHypothesis = hyp, alpha = sig)
      tobind <- data.frame(as.data.frame(thresh_res), 
                           type = hyp, 
                           comparison = paste(combos[[i]][1], combos[[i]][2], sep = "_vs_")) %>% 
        mutate(label = str_remove_all(comparison, "embryo|L1|L3")) %>%
        rownames_to_column(var = "WBGeneID")
      df <- rbind(df, tobind) 
    }
  }
  df <- df %>% 
    drop_na(padj) %>%
    mutate(isDE = case_when(padj < sig ~ TRUE,
                            padj > sig ~ FALSE,
                            is.na(padj) ~ FALSE))
  df
}


de_category_MA_plot <- function(df, title){
  df %>% filter(isDE == TRUE) %>%
    ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = type)) +
    geom_point(data =df %>% mutate(padj = replace_na(padj, 1)), shape = 16, alpha = 0.01, stroke = 0, size = 1, color = "grey") +
    geom_point(shape = 16, alpha = 0.5, stroke = 0, size = 1) +
    ylim(c(-10,10))+
    facet_wrap(~label) +
    # scale_color_manual(values = c("black", "red"), name = "q.value < 0.1") +
    theme_classic() +
    ggtitle(title)
}


de_category_bar_plot <- function(df, title){
  df %>% filter(isDE == TRUE) %>% group_by(label, type) %>% summarize(genes = n()) %>%
    ggplot(aes(x = type, y = genes, label = genes, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(vjust = -0.25) +
    facet_wrap(~label) +
    theme_classic() +
    ggtitle(title)
}



tissue_annotated_MA <- function(in_res, de_df){
  df <- as.data.frame(in_res) %>% rownames_to_column(var = "WBGeneID") %>%
    left_join(tissue_specific_genes, by = "WBGeneID") %>%
    mutate(tissue = replace_na(tissue, "other")) %>%
    left_join(de_df %>% filter(label == "GFPplus_vs_GFPminus") %>% select(WBGeneID, type, isDE), by = "WBGeneID")
  
  df %>% filter(isDE == TRUE) %>%
    ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = type)) +
    geom_point(data =df %>% select(-tissue), shape = 16, alpha = 0.1, stroke = 0, size = 1, color = "grey") +
    geom_point(shape = 16, alpha = 0.5, stroke = 0, size = 1) +
    facet_wrap(~tissue) +
    ylim(c(-10,10)) +
    theme_classic()
}

tissue_gene_quant <- function(in_df, sig = 0.01, thresh = 1){
  my_plot <- in_df %>% filter(isDE == TRUE, label == "GFPplus_vs_GFPminus") %>%
    left_join(tissue_specific_genes, by = "WBGeneID") %>%
    mutate(tissue = replace_na(tissue, "other"), padj = replace_na(padj, 1)) %>%
    group_by(tissue, type) %>%
    summarise(genes = n()) %>%
    ggplot(aes(x = type, y = genes, label = genes, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(vjust = -0.25) +
    facet_wrap(~tissue)+
    ggtitle(paste("comparison: ",deparse(substitute(in_df)), "\nlfc = ",thresh," & padj < ",sig, sep = "")) +
    theme_classic()
  my_plot
}