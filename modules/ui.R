### DEAKit - User Interface
### Made by: Pablo Rodríguez

# Load modules
for(m in modules.tabs){
  print(paste("Loading ui module: ", m, sep=''))
  source(paste("modules/", m, "/", m, "Panel.R", sep=''), local = TRUE)
}

# Define UI for application that draws a histogram
app.ui <- shinyUI(fluidPage(
  fluidPage(title = "DEAKit - Application",
            theme = "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css",
            tags$script(src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"),
            tags$script("$(document).ready(function() { $(\"body\").tooltip({ selector: '[data-toggle=tooltip]' });});"),
            useShinyjs(),
    
            # Page
            navbarPage(
              div(
                img(
                  src = "https://vignette.wikia.nocookie.net/fantendo/images/5/5c/Corruption_and_Energy.png/revision/latest?cb=20151117140508",#"R/www/logo.png",
                  height = 30,
                  width = "auto"
                  )
                ),
              
              load.data.panel,
              quality.control.panel,
              normalization.panel,
              dea.panel,
              pathview.panel,
              #),
            
              includeScript("www/js/mainFunctions.js")
            )
  )
  )
)
