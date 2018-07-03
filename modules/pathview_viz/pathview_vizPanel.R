pathview.panel <- tabPanel(title = "Pathway Enrichment",
                      fluidRow(
                        column(12,
                               span("Pathway enrichemt analysis parameters:"),
                               wellPanel(
                                 div(class = "col-md-12",
                                     #div(class = "col-md-6",
                                     textInput("pathview.all.path", "Pathway enrichemt analysis output of all the significantly expressed genes:", value = file.path(getwd(), base.output.path, "PCR-0-1_Enrichment-ALL.csv"),
                                               placeholder = "/home/userXY/output/normalized-metadata.rds"),
                                     helpText("Path of the output for the pathway enrichment analysis of all significantly expressed genes"),
                                     
                                     textInput("pathview.over.path", "Pathway enrichemt analysis output of the significantly over-expressed genes:", value = file.path(getwd(), base.output.path, "PCR-0-1_Enrichment-OVER.csv"),
                                               placeholder = "/home/userXY/output/normalized-metadata.rds"),
                                     helpText("Path of the output for the pathway enrichment analysis of significantly over-expressed genes"),
                                     
                                     textInput("pathview.under.path", "Pathway enrichemt analysis output of the significantly under-expressed genes:", value = file.path(getwd(), base.output.path, "PCR-0-1_Enrichment-UNDER.csv"),
                                               placeholder = "/home/userXY/output/normalized-metadata.rds"),
                                     helpText("Path of the output for the pathway enrichment analysis of significantly under-expressed genes"),
                                     
                                     # LogFC and Pvalue Thresholds
                                     #div(class = "col-md-6", 
                                     numericInput("p.value.pathview", "P Value", value = 0.05, min=0, max=1),
                                     helpText("P Value threshold used for pathway enrichment analysis")
                                 ),
                                 
                                 div(class="text-center", 
                                     useShinyjs(),
                                     tags$style(appCSS),
                                     withBusyIndicatorUI(
                                       actionButton("pathviewButton", 
                                                    "Execute pathway enrichment analysis", 
                                                    class="btn-primary"))
                                 )
                               )
                        ),
                        column(12,
                                span("Enriched pathways dotplots:"),
                                wellPanel(id = "enr.all.dot",
                                          div(class = "col-md-12 center", 
                                            span("Dotplot of enriched pathways from all significantly expressed genes:"),
                                            plotOutput("all.enr.dotplot"),
                                 
                                            span("Dotplot of enriched pathways from the significantly over-expressed genes:"),
                                            plotOutput("over.enr.dotplot"),
                                            
                                            span("Dotplot of enriched pathways from the significantly under-expressed genes:"),
                                            plotOutput("under.enr.dotplot")
                                          
                                )
                        )),
                        column(12,
                               span("Pathway filtering:"),
                               wellPanel(
                                 div(class = "col-md-12",
                                     #div(class = "col-md-6",
                                     numericInput("q.value.pathview", "Q-value", value = 1e-10, min=0, max=1),
                                     helpText("Q-value threshold for filtering the significantly enriched pathways"),
                                     
                                     selectInput("filter.pathview", "Select a subset of significant genes:", choices = c("All", "Over", "Under"),
                                                 selected = "All"),
                                     helpText("Indicate if the pathways are filtered from the enrichment of all, the significantly over or under expressed genes"),
                                     
                                     HTML("<b>Select if pathways are also filtered from an input file:</b>"),
                                     checkboxInput("check_filter_file", "Filter pathways from text file"),
                                     
                                     conditionalPanel(
                                       condition = "input.check_filter_file == true",
                                       textInput("pathview.filter.file", "Text file containing the filtered pathways:", value = file.path(getwd(), "data", "PanCancer_Pathways.txt")),
                                       helpText("Path of the input file containing the Kegg IDs of the filtered pathways, one for each line")
                                     )
                                ),
                                
                                div(class="text-center", 
                                    useShinyjs(),
                                    tags$style(appCSS),
                                    withBusyIndicatorUI(
                                      actionButton("pathviewQvalueButton", 
                                                   "Execute q-value pathway filtering", 
                                                   class="btn-primary"))
                                )
                               )
                        ),
                        column(12,
                               span("Filtered enriched pathways:"),
                               wellPanel(id = "path.res",
                                         div(class = "col-md-12",
                                             span("Select one of the pathways:"),
                                             dataTableOutput("loaded.enrich")
                              
                                        ),
                                        
                                        div(class="text-center", 
                                            useShinyjs(),
                                            tags$style(appCSS),
                                            withBusyIndicatorUI(
                                              actionButton("pathviewImageButton", 
                                                           "View Kegg pathway", 
                                                           class="btn-primary"))
                                        )
                               )
                               
                        ),
                        fillCol(width = "100%", height = "100%", #column(12,
                               #span("View Kegg selected pathway:"),
                               #wellPanel(id = "path.res",
                                         div(class = "col-md-12",
                                             span("Selected pathway image:"),
                                             imageOutput("kegg.image")#, width = PLOT.SIZE, height = PLOT.SIZE)
                                             
                                         )
                               #)
                               
                        )
                        
                )

)