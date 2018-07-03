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
source('DEA_functions.R')
source('helpers.R')

# Setup
SEED=12345
set.seed(SEED)

# Output directories, only the name of the output files, no final '/'
base.output.path <- "output-app"
if(!file.exists(base.output.path)) dir.create(base.output.path, showWarnings = FALSE)

img.output.path <- NULL
if(is.null(img.output.path)) {
  img.output.path <- file.path(base.output.path, "images")
  dir.create(img.output.path, showWarnings = FALSE)
} else if(!file.exists(img.output.path)) dir.create(img.output.path)


# Load modules
source('modules/modules.R', local = TRUE)
source('modules/ui.R', local = TRUE)
source('modules/server.R', local = TRUE)

# Run app
shinyApp(
  ui = app.ui,
  server = app.server
)

