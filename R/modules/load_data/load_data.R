load.data.values <- reactiveValues()

load.data.core <- eventReactive(input$loadDataButton, {
  rcc.path <- input$rcc.path
  rlf.path <- input$rlf.path   
  key.of.interest <- input$key.of.interest
  label.of.key <- input$label.of.key
  key.value.1 <- input$key.value.1
  key.value.2 <- input$key.value.2
  
  # Setup
  
  load.data.values$keys.vector <- unlist(lapply(seq_len(input$key.of.interest), function(x) paste("key",x,sep = "")), use.names = FALSE)
  load.data.values$values.of.key.of.interest <- c(key.value.1, key.value.2)
  
  # Load data + Extract metadata
  #output$metadata <- extract.rcc.metadata(rcc.path, keys.vector, values.of.key.of.interest)
  load.data.values$metadata <- extract.rcc.metadata(rcc.directory, keys.vector, VALUES.OF.KEY.OF.INTEREST)
  
  # Retrieve Set and Count matrix
  #rcc.set.and.count.matrix <- get.RCC.set.and.counts(load.data.values$metadata, rcc.path, rlf.path)
  load.data.values$rcc.set.and.count.matrix <- get.RCC.set.and.counts(metadata, rcc.directory, rlf.filename)
  
  load.data.values$eset <- rcc.set.and.count.matrix$set
  
  load.data.values$counts <- rcc.set.and.count.matrix$count.matrix
  
  #shinyjs::enable(id = "searchButton")
  
  return(as.data.frame(load.data.values$metadata))
})

output$loaded.metadata <- renderDataTable({
  return(load.data.core())
}, options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE), escape=FALSE, selection = 'single'
)

output$loaded.counts <- renderDataTable({
  return(load.data.values$counts[,c(1,2,3,4,5)])
}, options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE), escape=FALSE, selection = 'single'
)