pathview.values <- reactiveValues()

pathview.core <- function(){
  pathview.values$p.value <- input$p.value.pathview
  shiny.PATHVIEW.ALL.OUTPUT.NAME <- input$pathview.all.path
  shiny.PATHVIEW.OVER.OUTPUT.NAME <- input$pathview.over.path
  shiny.PATHVIEW.UNDER.OUTPUT.NAME <- input$pathview.under.path
  
  # The core script for the pathway analysis expects a dataframe with the DEA named res
  res <- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) 
  
  # Setup ENTREZ - Takes time
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
  entrez <- biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
  entrez$entrezgene <- as.character(entrez$entrezgene)
  entrezsymbol <- biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
  entrezsymbol$entrezgene <- as.character(entrezsymbol$entrezgene)
  
  converted <- unlist(lapply(strsplit(res$rowname, "_", fixed = TRUE), "[", 4))
  converted <- unlist(lapply(strsplit(converted, ".", fixed = TRUE), "[", 1))
  converted <- paste0("NM_", converted)
  res$refseq_mrna <- converted
  
  # Enrich - All
  enrich.rs <- enrich.cp(res, entrez, label.of.interest, type="all", pval.threshold = pathview.values$p.value, lfc.threshold = log.fc)
  enrich.rs.summary <- enrich.rs$summary %>% arrange(p.adjust)
  enrich.rs.summary <- convert.enriched.ids(enrich.rs.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
  
  # Enrich - Over
  enrich.rs.over <- enrich.cp(res, entrez, label.of.interest, type="over", pval.threshold = pathview.values$p.value, lfc.threshold = log.fc)
  enrich.rs.over.summary <- enrich.rs.over$summary %>% arrange(p.adjust)
  enrich.rs.over.summary <- convert.enriched.ids(enrich.rs.over.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
  
  # Enrich - Under
  enrich.rs.under <- enrich.cp(res, entrez, label.of.interest, type="under", pval.threshold = pathview.values$p.value, lfc.threshold = log.fc)
  enrich.rs.under.summary <- enrich.rs.under$summary %>% arrange(p.adjust)
  enrich.rs.under.summary <- convert.enriched.ids(enrich.rs.under.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
  
  write_csv(x = enrich.rs.summary, path = shiny.PATHVIEW.ALL.OUTPUT.NAME)
  write_csv(x = enrich.rs.over.summary, path = shiny.PATHVIEW.OVER.OUTPUT.NAME)
  write_csv(x = enrich.rs.under.summary, path = shiny.PATHVIEW.UNDER.OUTPUT.NAME)
  
  
  ### Dotplots
  DOTPLOT.ENRICHMENT.ALL.TITLE <- "Dotplot of All Enriched Pathways after DEA"
  DOTPLOT.ENRICHMENT.OVER.TITLE <- "Dotplot of Overexpressed Enriched, Pathways after DEA"
  DOTPLOT.ENRICHMENT.UNDER.TITLE <- "Dotplot of Underexpressed Enriched Pathways after DEA"
  
  # KEGG enrichment in all significantly expressed genes:
  pathview.values$all.dotplot <- dotplot(enrich.rs$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.ALL.TITLE)
  # KEGG enrichment in over-expressed genes:
  pathview.values$over.dotplot <- dotplot(enrich.rs.over$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.OVER.TITLE)
  # KEGG enrichment in under-expressed genes
  pathview.values$under.dotplot <- dotplot(enrich.rs.under$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.UNDER.TITLE)
  
  
  # OUTPUT
  output$all.enr.dotplot <- renderPlot(plot(pathview.values$all.dotplot))
  output$over.enr.dotplot <- renderPlot(plot(pathview.values$over.dotplot))
  output$under.enr.dotplot <- renderPlot(plot(pathview.values$under.dotplot))
  
  # Store variables used in other scripts in Global Environment
  assign("enrich.rs.summary", enrich.rs.summary, envir = globalenv())
  assign("enrich.rs.over.summary", enrich.rs.over.summary, envir = globalenv())
  assign("enrich.rs.under.summary", enrich.rs.under.summary, envir = globalenv())
  
}

pathway.qvalue.core <- function() {
  pathview.values$q.value <- input$q.value.pathview
  enrich.summary <- switch(input$filter.pathview,
                           All = enrich.rs.summary,
                           Over = enrich.rs.over.summary,
                           Under = enrich.rs.under.summary)
  
  enrich.kegg.summary.filtered <- ontology.enrichment.qval.filter(enrichment.summary = enrich.summary,
                                                                      ontology = "kg",
                                                                      qval = pathview.values$q.value)
  assign("enrich.kegg.summary.filtered", enrich.kegg.summary.filtered, envir = globalenv())
  
  # OUTPUT
  output$loaded.enrich <- renderDataTable(enrich.kegg.summary.filtered[, 1:7], 
                                            options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE, scrollX = TRUE), escape=FALSE, selection = 'single'
  )
}

pathway.view.core <- function() {
  pathview.values$selected.row <- input$loaded.enrich_rows_selected
  pathview.values$selected.pathway <- enrich.kegg.summary.filtered[pathview.values$selected.row, ]$ID
  
  
  # Create PathviewR pathway
  setwd("images")
  viewkeggpath(path=pathview.values$selected.pathway, enrichment = enrich.rs.summary, dea.p = dea.p, output = ".")
  setwd("..")
  
  # Output
  output$kegg.image <- renderImage(list(src = paste0("images/", pathview.values$selected.pathway, ".pathview.png")))
}


observeEvent(input$pathviewButton, pathview.core())

observeEvent(input$pathviewQvalueButton, pathway.qvalue.core())

observeEvent(input$pathviewImageButton, pathway.view.core())

