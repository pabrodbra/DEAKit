### ------------------------
### pathfind.R
### Author: Pablo Rodriguez Brazzarola & Guillermo Lopez Garcia
### Description: Script to manually perform the pathfindR Pathway Enrichment Analysis, 
###              described in https://github.com/egeulgen/pathfindR
### ------------------------

### Load functions
source("DEA_functions.R")

### Test using default data
data("RA_input")
RA_input

# Execute Enrichment Workflow
setwd("pathfindR")
test.result <- run_pathfindR(RA_input, p_val_threshold = 0.05, enrichment_threshold = 0.05,
                        adj_method = "BH", search_method = "GR", iterations = 1, 
                        pin_name_path = "Biogrid")

# Execute Clustering Workflow  
choose_clusters(test.result)
setwd("../")


### Using our own data
# Input data needs to be in the format produced by our own getDE.raw function
# We need to generate this csv for the particular DEA analysis
dea.raw <- read.csv("pathfindR/PCR_DEA_raw.csv", stringsAsFactors = FALSE)

# Pre-process the data to be in the correct format 
gds <- gene.dea.summary(dea.raw, adj.pval.threshold = 0.05)

# Execute Enrichment Workflow
setwd("pathfindR")
result <- run_pathfindR(gds, p_val_threshold = 0.05, enrichment_threshold = 0.05,
                        adj_method = "BH", search_method = "GR", iterations = 5, 
                        pin_name_path = "Biogrid")

# Execute Clustering Workflow (it creates a shiny HTML document)
choose_clusters(result)
setwd("../")
