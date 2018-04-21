#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("R/DEA_functions.R")

# Load modules
source('R/modules/ui.R', local = TRUE)
source('R/modules/server.R', local = TRUE)

# Run app
shinyApp(
  ui = app.ui,
  server = app.server
)