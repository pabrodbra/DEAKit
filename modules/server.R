### DEAKit - Server
### Made by: Pablo Rodríguez

app.server <- shinyServer(function(input, output, session) {
  
  # Load modules
  for(m in modules.tabs){
    print(paste("Loading server module: ", m, sep=''))
    source(paste("modules/", m, "/", m, ".R", sep=''), local = TRUE)
  }
  
})
