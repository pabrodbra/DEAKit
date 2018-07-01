load.data.values <- reactiveValues()

load.data.core <- function() {
  rcc.path <- input$rcc.path
  rlf.path <- input$rlf.path   
  key.of.interest <- input$key.of.interest
  label.of.key <- input$label.of.key
  key.value.1 <- input$key.value.1
  key.value.2 <- input$key.value.2
  
  # Setup
  load.data.values$keys.vector <- reactive(unlist(lapply(seq_len(key.of.interest), function(x) paste("key",x,sep = "")), use.names = FALSE))
  load.data.values$key.label <- reactive(tail(load.data.values$keys.vector(), 1))
  load.data.values$values.of.key.of.interest <- reactive(c(key.value.1, key.value.2))
  load.data.values$label.of.interest <- reactive(label.of.key)
  
  # Load data + Extract metadata
  load.data.values$metadata <- reactive(extract.rcc.metadata(rcc.path, load.data.values$keys.vector(), load.data.values$values.of.key.of.interest()))
  
  # Retrieve Set and Count matrix
  load.data.values$rcc.set.and.count.matrix <- get.RCC.set.and.counts(load.data.values$metadata(), rcc.path, rlf.path)
  load.data.values$eset <- reactive(load.data.values$rcc.set.and.count.matrix$set)
  load.data.values$counts <- reactive(load.data.values$rcc.set.and.count.matrix$count.matrix)
  
  # Store variables used in other scripts in Global Environment
  assign("key.label", load.data.values$key.label(), envir = globalenv())
  assign("label.of.interest", load.data.values$label.of.interest(), envir = globalenv())
  assign("metadata", load.data.values$metadata(), envir = globalenv())
  assign("eset", load.data.values$eset(), envir = globalenv())
  assign("counts", load.data.values$counts(), envir = globalenv())
  
  # OUTPUT
  output$loaded.metadata <- renderDataTable(load.data.values$metadata(), 
    options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE), escape=FALSE, selection = 'single'
  )
  
  output$loaded.counts <- renderDataTable(load.data.values$counts()[,c(1,2,3,4,5)],
    options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE), escape=FALSE, selection = 'single'
  )
}

observeEvent(input$loadDataButton, load.data.core())

#hideElement("ld.res")