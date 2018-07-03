normalization.values <- reactiveValues()

normalization.core <- function() {
  withBusyIndicatorServer("normalizeButton", {
    shiny.RDS.METADATA.OUTPUT <- input$norm.metadata.path
    shiny.RDS.NORMALIZED.OUTPUT <- input$norm.count.path
    PRINCIPAL.COMPONENT.X <- input$pca.x
    PRINCIPAL.COMPONENT.Y <- input$pca.y
    
    FOV.THRESHOLD <- input$fov.threshold
    LOD.THRESHOLD <- (100-FOV.THRESHOLD)/100 # Dependant of FOV.THRESHOLD
    
    # Pre-Normalization
    metadata$pos_nf <- pos.factor(eset) # Positive Control normalization
    counts <- counts[!grepl("Positive", rownames(counts)),]
    counts <- counts[!grepl("Negative", rownames(counts)),]
    
    # Plot pre-normalization
    pre.normalization.boxplot <- boxplot.count.matrix(counts, "Pre-Normalization with Positive controls")
    normalization.values$pre.norm.plot <- pre.normalization.boxplot + geom_hline(yintercept = lod, colour="red")
    
    # Normalize
    normalization.values$ncounts <- counts %*% diag(metadata$pos_nf) # Positive Control normalization
    colnames(normalization.values$ncounts) <- colnames(counts)
    
    # Plot post-normalization
    post.normalization.boxplot <- boxplot.count.matrix(normalization.values$ncounts, "Post-Normalization with Positive controls")
    normalization.values$post.norm.plot <- post.normalization.boxplot + geom_hline(yintercept = lod, colour="red")
    # ^^ No noticeable change ^^ #
    
    # Housekeeping normalization
    # Normalize housekeeping
    metadata$hk_nf <- hk.factor(counts, lod) # Housekeeping normalization
    normalization.values$ncounts <- counts %*% diag(metadata$hk_nf) # Housekeeping normalization
    colnames(normalization.values$ncounts) <- colnames(counts)
    
    # Plot normalized housekeeping
    housekeeping.post.normalization.boxplot <- boxplot.count.matrix(normalization.values$ncounts, "Post-Normalization with Housekeeping controls")
    normalization.values$norm.plot <- housekeeping.post.normalization.boxplot + geom_hline(yintercept = lod, colour="red") + geom_smooth(se=T, aes(group=1))
    
    # Drop genes below threshold
    all.names <- rownames(normalization.values$ncounts)
    normalization.values$ncounts <- normalization.values$ncounts[(rowSums(normalization.values$ncounts < lod) < round((LOD.THRESHOLD * ncol(normalization.values$ncounts)),0)),]
    filtered.names <- rownames(normalization.values$ncounts)
    normalization.values$filter.out.genes <- all.names[all.names %!in% filtered.names]
    
    # Save metadata dataframe and the normalised counts matrix - CHANGE FOR INPUT
    saveRDS(metadata, shiny.RDS.METADATA.OUTPUT) 
    saveRDS(normalization.values$ncounts, shiny.RDS.NORMALIZED.OUTPUT)
    
    # PCA
    #PCA varianza
    pcdata <- scale(normalization.values$ncounts, center = TRUE, scale = TRUE)
    pc <- pca.loadings(pcdata, 50)
    comps <- data.frame(pc$x)
    comps$Name <- rownames(comps)
    comps <- comps %>% left_join(metadata, by = c(Name = "Sample.name"))
    
    # Pass labels to legend and/or title
    
    PCA.PLOT.TITLE <- paste("PCA Plot | Component", PRINCIPAL.COMPONENT.X, "-", PRINCIPAL.COMPONENT.Y, sep = " ")
    
    normalization.values$pca.plot <- plotPCA(comps, pc, PRINCIPAL.COMPONENT.X, PRINCIPAL.COMPONENT.Y, key.label, 
                          legend.label =  label.of.interest, plot.title = PCA.PLOT.TITLE)
    
    # Store variables used in other scripts in Global Environment
    assign("metadata", metadata, envir = globalenv())
    assign("counts", counts, envir = globalenv())
    assign("ncounts", normalization.values$ncounts, envir = globalenv())
    
    # OUTPUT
    output$pre.norm.plot <- renderPlot(plot(normalization.values$pre.norm.plot))
    pdf(file.path(base.output.path, "Pre-Normalization.pdf"))
    plot(normalization.values$pre.norm.plot)
    invisible(dev.off())
    
    output$post.norm.plot <- renderPlot(plot(normalization.values$post.norm.plot))
    pdf(file.path(base.output.path, "Positive-Normalization.pdf"))
    plot(normalization.values$post.norm.plot)
    invisible(dev.off())
    
    output$norm.plot <- renderPlot(plot(normalization.values$norm.plot))
    pdf(file.path(base.output.path, "Housekeeping-Normalization.pdf"))
    plot(normalization.values$norm.plot)
    invisible(dev.off())
    
    output$pca.plot <- renderPlot(normalization.values$pca.plot)
    pdf(file.path(base.output.path, "Samples-PCA.pdf"))
    plot(normalization.values$pca.plot)
    invisible(dev.off())
  })
}

observeEvent(input$normalizeButton, normalization.core())
