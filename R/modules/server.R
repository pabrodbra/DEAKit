### DEAKit - Server
### Made by: Pablo Rodríguez

app.server <- shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
  # Load modules
  for(m in modules.tabs){
    print(paste("Loading server module: ", m, sep=''))
    source(paste("R/modules/", m, "/", m, ".R", sep=''), local = TRUE)
  }
  
})
