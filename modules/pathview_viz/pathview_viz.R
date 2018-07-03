pathview.values <- reactiveValues()

pathview.core <- function() {
  withBusyIndicatorServer("pathviewButton", {
    pathview.values$p.value <- input$p.value.pathview
    shiny.PATHVIEW.ALL.OUTPUT.NAME <- input$pathview.all.path
    shiny.PATHVIEW.OVER.OUTPUT.NAME <- input$pathview.over.path
    shiny.PATHVIEW.UNDER.OUTPUT.NAME <- input$pathview.under.path
    
    # The core script for the pathway analysis expects a dataframe with the DEA named res
    res <- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) 
    
    # Setup ENTREZ - Takes time
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = biomart_dataset)
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
    pdf(file.path(base.output.path, "Dotplot-All.pdf"))
    plot(pathview.values$all.dotplot)
    invisible(dev.off())
    
    output$over.enr.dotplot <- renderPlot(plot(pathview.values$over.dotplot))
    pdf(file.path(base.output.path, "Dotplot-Over.pdf"))
    plot(pathview.values$over.dotplot)
    invisible(dev.off())
    
    output$under.enr.dotplot <- renderPlot(plot(pathview.values$under.dotplot))
    pdf(file.path(base.output.path, "Dotplot-Under.pdf"))
    plot(pathview.values$under.dotplot)
    invisible(dev.off())
    
    # Store variables used in other scripts in Global Environment
    assign("enrich.rs.summary", enrich.rs.summary, envir = globalenv())
    assign("enrich.rs.over.summary", enrich.rs.over.summary, envir = globalenv())
    assign("enrich.rs.under.summary", enrich.rs.under.summary, envir = globalenv())
  
  })
}

pathway.qvalue.core <- function() {
  withBusyIndicatorServer("pathviewQvalueButton", {
    pathview.values$q.value <- input$q.value.pathview
    enrich.summary <- switch(input$filter.pathview,
                             All = enrich.rs.summary,
                             Over = enrich.rs.over.summary,
                             Under = enrich.rs.under.summary)
    
    # Filter using Q-value
    enrich.kegg.summary.filtered <- ontology.enrichment.qval.filter(enrichment.summary = enrich.summary,
                                                                        ontology = "kg",
                                                                        qval = pathview.values$q.value)
    # Filter using input file
    if(input$check_filter_file) {
      filter.path <- scan(input$pathview.filter.file, what = typeof(enrich.kegg.summary.filtered$ID), 
                          sep = "\n", quiet = TRUE)
      enrich.kegg.summary.filtered <- enrich.kegg.summary.filtered %>% filter(ID %in% filter.path)
    }
    
    assign("enrich.kegg.summary.filtered", enrich.kegg.summary.filtered, envir = globalenv())
    
    # OUTPUT
    output$loaded.enrich <- renderDataTable(enrich.kegg.summary.filtered[, 1:7], 
                                              options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE, scrollX = TRUE), escape=FALSE, selection = 'single'
    )
  })
}

pathway.view.core <- function() {
  withBusyIndicatorServer("pathviewImageButton", {
    pathview.values$selected.row <- input$loaded.enrich_rows_selected
    
    # Check a row is selected
    if(!length(pathview.values$selected.row))
      stop("You must click on a Pathway first!")
    
    # Extract selected pathway ID
    pathview.values$selected.pathway <- enrich.kegg.summary.filtered[pathview.values$selected.row, ]$ID
    
    # Create PathviewR pathway
    curr.dir <- getwd()
    message(curr.dir)
    setwd(img.output.path)
    viewkeggpath(path=pathview.values$selected.pathway, enrichment = enrich.kegg.summary.filtered, dea.p = dea.p, output = ".")
    setwd(curr.dir)
    
    # Output
    output$kegg.image <- renderImage(list(src = file.path(img.output.path, paste0(pathview.values$selected.pathway, ".pathview.png"))),
                                     deleteFile = FALSE)
  })
}

observeEvent(input$pathviewButton, pathview.core())

observeEvent(input$pathviewQvalueButton, pathway.qvalue.core())

observeEvent(input$pathviewImageButton, pathway.view.core())

