quality.control.values <- reactiveValues()

quality.control.core <- function(){
  fov.threshold <- input$fov.threshold
  bd.threshold.min <- input$bd.threshold.min
  bd.threshold.max <- input$bd.threshold.max
  bd.thresholds <- c(bd.threshold.min, bd.threshold.max)
  
  # Field of View (FOV) Plots
  quality.control.values$fov <- plotFOV(eset = eset, metadata = metadata, fov.threshold = fov.threshold,
                                        comparison.key = key.label, legend.label = label.of.interest)
  
  # Binding Density (BD) Plots
  quality.control.values$bd <- plotBD(eset = eset, metadata = metadata, y.thresholds = bd.thresholds, 
                                      comparison.key = key.label, legend.label = label.of.interest)
  
  # Positive Controls
  quality.control.values$positive.control <- boxplot.expr(eset, is.positive)
  
  # Negative Controls
  quality.control.values$negative.control <- boxplot.expr(eset, is.negative)
  
  # Noise Threshold
  quality.control.values$lodcounts <- extract.pred(eset, is.negative)
  quality.control.values$lod <- mean(quality.control.values$lodcounts$count) + 2 * sd(quality.control.values$lodcounts$count)
  
  # Housekeeping Genes
  quality.control.values$housekeeping.boxplot <- boxplot.expr(eset, is.housekeeping)
  quality.control.values$housekeeping.boxplot <- quality.control.values$housekeeping.boxplot+geom_hline(yintercept = (quality.control.values$lod),colour="red")
  
  # Expression of all the housekeeping genes in each sample
  quality.control.values$housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Housekeeping"), 
                                                                title.label = "Housekeeping Genes Expression per Sample")
  quality.control.values$housekeeping.condition.boxplot <- quality.control.values$housekeeping.condition.boxplot+geom_hline(yintercept = (quality.control.values$lod),colour="red")
  
  # Expression of multiple condition genes in each sample
  quality.control.values$endogenous.housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Endogenous", "Housekeeping"), 
                                                                           title.label = "Housekeeping + Endogenous Genes Expression per Sample")
  quality.control.values$endogenous.housekeeping.condition.boxplot <- quality.control.values$endogenous.housekeeping.condition.boxplot+geom_hline(yintercept = (quality.control.values$lod),colour="red")
  
  # Store variables used in other scripts in Global Environment
  assign("lodcounts", quality.control.values$lodcounts, envir = globalenv())
  assign("lod", quality.control.values$lod, envir = globalenv())
  
  #showElement("qc.res")
  # OUTPUT
  output$qc.fov <- renderPlot(plot(quality.control.values$fov))
  output$qc.bd <- renderPlot(plot(quality.control.values$bd))
  output$qc.positive.control.bp <- renderPlot(plot(quality.control.values$positive.control))
  output$qc.negative.control.bp <- renderPlot(plot(quality.control.values$negative.control))
  output$qc.housekeeping.bp <- renderPlot(plot(quality.control.values$housekeeping.boxplot))
  output$qc.sample.housekeeping.bp <- renderPlot(plot(quality.control.values$housekeeping.condition.boxplot))
  output$qc.sample.endogenous.housekeeping.bp <- renderPlot(plot(quality.control.values$endogenous.housekeeping.condition.boxplot))
  
}

observeEvent(input$qualityControlButton, quality.control.core())

#hideElement("qc.res")