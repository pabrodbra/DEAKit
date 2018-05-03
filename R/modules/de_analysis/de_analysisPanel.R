dea.panel <- tabPanel(title = "DEA",
            fluidRow(
              column(12,
                     span("Differential Expression Analysis parameters:"),
                     wellPanel(
                       div(class = "col-md-12",
                           #div(class = "col-md-6",
                           textInput("dea.path", "DEA Enrichment output:",
                                     placeholder = "/home/userXY/output/normalized-metadata.rds"),
                           helpText("Path of the output for the DEA enrichment data"),
                           
                           # RDS - Normalized count matrix output
                           textInput("dea.sign.path", "DEA Significant Enrichment output:",
                                     placeholder = "/home/userXY/output/normalized-count.rds"),
                           helpText("Path of the output for the DEA Significant Enrichment data")
                           ,
                           #),
                           
                           # LogFC and Pvalue Thresholds
                           #div(class = "col-md-6", 
                           numericInput("p.value.dea", "P Value", value = 0.05, min=0, max=1),
                           helpText("P Value threshold"),
                           
                           numericInput("log.fc.dea", "Log FC", value = 1, min=0, max=10),
                           helpText("Log Fold Change threshold")
                           
                       ),
                       
                       div(class="text-center", actionButton("deaButton", "Execute Differential Expression Analysis", class="btn-primary"))
                     )
              ),
              column(12, 
                     span("Differential Expression Analysis results:"),
                     wellPanel(id = "dea.res",
                               span("Volcano Plot:"),
                               plotOutput("volcano.plot",  width = PLOT.SIZE, height = PLOT.SIZE)
                     )
              )
            )
)