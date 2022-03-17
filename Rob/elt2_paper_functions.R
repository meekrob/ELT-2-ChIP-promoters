
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