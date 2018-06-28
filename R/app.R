#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(DT)
#('shiny', 'shinyjs', 'DT')
source('R/DEA_functions.R')

# Load modules
source('R/modules/modules.R', local = TRUE)
source('R/modules/ui.R', local = TRUE)
source('R/modules/server.R', local = TRUE)

# Setup
SEED=12345
set.seed(SEED)
base.output.path <- "output/"

# Run app
shinyApp( 
  ui = app.ui,
  server = app.server
)

