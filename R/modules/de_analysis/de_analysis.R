dea.values <- reactiveValues()

dea.core <- function(){
  dea.values$p.value <- input$p.value.dea
  dea.values$log.fc <- input$log.fc.dea
  shiny.DEA.OUTPUT.NAME <- input$dea.path
  shiny.DEA.SIGNIFICATIVE.OUTPUT.NAME <- input$dea.sign.path
  
  design <- as.formula(paste("~", key.label, sep = " "))
  rownames(normalization.values$metadata) <- normalization.values$metadata$Sample.name
  
  dea.values$dea.o <- getDE(ncounts = round(ncounts), metadata = normalization.values$metadata, design = design)
  dea.values$dea.p <- dea.o %>% dplyr::filter(padj < dea.values$p.value) %>% arrange(-log2FoldChange)
  
  # Save DEA Results
  
  write_csv(dea.values$dea.o, DEA.OUTPUT.NAME)
  write_csv(dea.values$dea.p, DEA.SIGNIFICATIVE.OUTPUT.NAME)
  
  # Show DEA Results
  #dplyr::select(dea.p, type, gene.name, id, baseMean, log2FoldChange, pvalue, padj)
  
  # Volcano Plot
  tab <- data.frame(logFC = dea.values$dea.o$log2FoldChange, negLogPval = -log10(dea.values$dea.o$padj), Gene=dea.values$dea.o$gene.name)
  
  dea.values$volcano.plot <- plot.Volcano2(tab, dea.values$log.fc, dea.values$p.value, plot.title = "Volcano Plot")
  
  # OUTPUT
  output$volcano.plot <- renderPlot(plot(dea.values$volcano.plot))
}

observeEvent(input$deaButton, dea.core())

