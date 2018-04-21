### DEAKit - User Interface
### Made by: Pablo Rodríguez

# Load modules
for(m in modules.tabs){
  print(paste("Loading ui module: ", m, sep=''))
  source(paste("R/modules/", m, "/", m, "Panel.R", sep=''), local = TRUE)
}

# Define UI for application that draws a histogram
app.ui <- shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 50,
                   value = 30)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot")
    )
  )
))
