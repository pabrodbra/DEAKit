### DEAKit - Modules helper
### Made by: Pablo Rodríguez
PLOT.SIZE <- "600px"
# Modules to be loaded.
# Each module must consists of one directory with two files: 
#     ./R/modules/<module_name>
#       |-- <module_name>.R
#       '-- <module_name>Panel.R.
modules.tabs <- c("load_data", "quality_control", "normalization" , "de_analysis", "pathview_viz")#, "pathfindR_viz")

# Modules to be disabled until search is perform.
# This list consists of the values of the tab panels.
#     e.g. tabPanel("Visualization", value = "visTab", ...)
#modules.hidden.tabs <- c("quality_control" , "de_analysis", "pathway_visualization")

# Proxy settings
#Sys.setenv(http_proxy="")