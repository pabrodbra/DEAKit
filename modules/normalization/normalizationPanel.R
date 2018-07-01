normalization.panel <- tabPanel(title = "Normalization",
                                fluidRow(
                                  column(12,
                                         span("Normalization parameters:"),
                                         wellPanel(
                                           div(class = "col-md-12",
                                               #div(class = "col-md-6",
                                               # RDS - Normalized metadata output
                                               textInput("norm.metadata.path", "Normalized metadata output:", value = file.path(getwd(), "output/PCR-0-1_BMI_metadata.rds"),
                                                         placeholder = "/home/userXY/output/normalized-metadata.rds"),
                                               helpText("Path of the output for the normalized metadata"),
                                               
                                               # RDS - Normalized count matrix output
                                               textInput("norm.count.path", "Normalized count matrix output:", value = file.path(getwd(), "output/PCR-0-1_BMI_HK-normalized-counts.rds"),
                                                         placeholder = "/home/userXY/output/normalized-count.rds"),
                                               helpText("Path of the output for the normalized count matrix")
                                               ,
                                               #),
                                               
                                               #div(class = "col-md-6", 
                                               numericInput("pca.x", "PCA Component X", value = 1, min=1, max=4),
                                               helpText("Binding Density minimum threshold"),
                                               
                                               numericInput("pca.y", "PCA Component Y", value = 2, min=1, max=4),
                                               helpText("Binding Density maximum threshold")
                                           ),
                                           
                                           div(class="text-center", actionButton("normalizeButton", "Normalize metadata and count matrix", class="btn-primary"))
                                         )
                                  ),
                                  fillCol(width = "100%", height = "100%", # 12, 
                                         span("Normalization results results:"),
                                         wellPanel(id = "norm.res",
                                                   span("Pre-normalization with Positive Control genes:"),
                                                   plotOutput("pre.norm.plot"),
                                                   
                                                   span("Post-normalization with Positive Control genes:"),
                                                   plotOutput("post.norm.plot"),
                                                   
                                                   span("Post-normalization with Housekeeping genes:"),
                                                   plotOutput("norm.plot"),
                                                   
                                                   span("PCA Plot:"),
                                                   plotOutput("pca.plot")
                                         )
                                  )
                                )
)