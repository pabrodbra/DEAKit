normalization.values <- reactiveValues()

normalization.core <- function(){
  shiny.RDS.METADATA.OUTPUT <- input$norm.metadata.path
  shiny.RDS.NORMALIZED.OUTPUT <- input$norm.count.path
  PRINCIPAL.COMPONENT.X <- input$pca.x
  PRINCIPAL.COMPONENT.Y <- input$pca.y
  
  normalization.values$metadata <- load.data.values$metadata
  normalization.values$eset <- load.data.values$eset
  normalization.values$counts <- load.data.values$counts
  normalization.values$lodcounts <-  quality.control.values$lodcounts
  normalization.values$lod <- quality.control.values$lod 
  
  FOV.THRESHOLD <- input$fov.threshold
  LOD.THRESHOLD <- (100-FOV.THRESHOLD)/100 # Dependant of FOV.THRESHOLD
  
  # Pre-Normalization
  normalization.values$metadata$pos_nf = pos.factor(normalization.values$eset) # Positive Control normalization
  normalization.values$counts <- counts[!grepl("Positive", rownames(counts)),]
  normalization.values$counts <- counts[!grepl("Negative", rownames(counts)),]
  
  # Plot pre-normalization
  pre.normalization.boxplot <- boxplot.count.matrix(counts, "Pre-Normalization with Positive controls")
  normalization.values$pre.norm.plot <- pre.normalization.boxplot + geom_hline(yintercept = normalization.values$lod,colour="red")
  
  # Normalize
  normalization.values$ncounts <- normalization.values$counts %*% diag(normalization.values$metadata$pos_nf) # Positive Control normalization
  colnames(normalization.values$ncounts) <- colnames(normalization.values$counts)
  
  # Plot post-normalization
  post.normalization.boxplot <- boxplot.count.matrix(normalization.values$ncounts, "Post-Normalization with Positive controls")
  normalization.values$post.norm.plot <- post.normalization.boxplot + geom_hline(yintercept = normalization.values$lod,colour="red")
  # ^^ No noticeable change ^^ #
  
  # Housekeeping normalization
  # Normalize housekeeping
  normalization.values$metadata$hk_nf <- hk.factor(normalization.values$ncounts, normalization.values$lod) # Housekeeping normalization
  normalization.values$ncounts <- normalization.values$ncounts %*% diag(normalization.values$metadata$hk_nf) # Housekeeping normalization
  colnames(normalization.values$ncounts) <- colnames(normalization.values$counts)
  
  # Plot normalized housekeeping
  housekeeping.post.normalization.boxplot <- boxplot.count.matrix(normalization.values$ncounts, "Post-Normalization with Housekeeping controls")
  normalization.values$norm.plot <- housekeeping.post.normalization.boxplot + geom_hline(yintercept = normalization.values$lod,colour="red") + geom_smooth(se=T, aes(group=1))
  
  # Drop genes below threshold
  all.names <- rownames(ncounts)
  normalization.values$ncounts <- ncounts[(rowSums(ncounts < lod) < round((LOD.THRESHOLD * ncol(ncounts)),0)),]
  filtered.names <- rownames(ncounts)
  normalization.values$filter.out.genes <- all.names[all.names %!in% filtered.names]
  
  # Save metadata dataframe and the normalised counts matrix - CHANGE FOR INPUT
  saveRDS(normalization.values$metadata, RDS.METADATA.OUTPUT) 
  saveRDS(normalization.values$ncounts, RDS.NORMALIZED.OUTPUT)
  
  # PCA
  #PCA varianza
  pcdata <- scale(ncounts, center = TRUE, scale = TRUE)
  pc <- pca.loadings(pcdata, 50)
  comps <- data.frame(pc$x)
  comps$Name <- rownames(comps)
  comps <- comps %>% left_join(metadata, by = c(Name = "Sample.name"))
  
  #- Pass labels to legend and/or title
  
  PCA.PLOT.TITLE <- paste("PCA Plot | Component", PRINCIPAL.COMPONENT.X, "-", PRINCIPAL.COMPONENT.Y, sep = " ")
  
  normalization.values$pca.plot <- plotPCA(comps, PRINCIPAL.COMPONENT.X, PRINCIPAL.COMPONENT.Y, key.label
          , legend.label =  LABEL.OF.INTEREST, plot.title = PCA.PLOT.TITLE)
  
  # OUTPUT
  output$pre.norm.plot <- renderPlot(plot(normalization.values$pre.norm.plot))
  output$post.norm.plot <- renderPlot(plot(normalization.values$post.norm.plot))
  output$norm.plot <- renderPlot(plot(normalization.values$norm.plot))
  output$pca.plot <- renderPlot(plot(normalization.values$pca.plot))
}

observeEvent(input$normalizeButton, normalization.core())
