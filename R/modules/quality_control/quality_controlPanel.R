quality.control.panel <- tabPanel(title = "Quality Control",
                                  fluidRow(
                                    column(12,
                                           span("Quality Control parameters:"),
                                           wellPanel(
                                             div(class = "col-md-12",
                                                 #div(class = "col-md-6",
                                                 # 6 - FIELD OF VIEW THRESHOLD - For FOV plot
                                                 numericInput("fov.threshold", "Field of View threshold", value = 80, min=0, max=100),
                                                 helpText("Field of View threshold")
                                                 ,
                                                 #),
                                                 
                                                 # 7 - BINDING DENSITY THRESHOLDS (min and max) - For BD plot
                                                 #div(class = "col-md-6", 
                                                 numericInput("bd.threshold.min", "Binding Density minimum", value = 0.05, min=0, max=3),
                                                 helpText("Binding Density minimum threshold"),
                                                 
                                                 numericInput("bd.threshold.max", "Binding Density maximum", value = 2.25, min=0, max=3),
                                                 helpText("Binding Density maximum threshold")
                                                 
                                             ),
                                             
                                             div(class="text-center", actionButton("qualityControlButton", "Generate plots", class="btn-primary"))
                                           )
                                    ),
                                    column(12, 
                                           span("Quality Control results:"),
                                           wellPanel(id = "qc.res",
                                             span("Field of View:"),
                                             plotOutput("qc.fov",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Binding Density:"),
                                             plotOutput("qc.bd",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Positive Controls Boxplot:"),
                                             plotOutput("qc.positive.control.bp",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Negative Controls Boxplot:"),
                                             plotOutput("qc.negative.control.bp",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Housekeeping Genes Boxplot:"),
                                             plotOutput("qc.housekeeping.bp",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Housekeeping Genes per Sample Boxplot:"),
                                             plotOutput("qc.sample.housekeeping.bp",  width = PLOT.SIZE, height = PLOT.SIZE),
                                             
                                             span("Endogenous and Housekeeping Genes per Sample Boxplot:"),
                                             plotOutput("qc.sample.endogenous.housekeeping.bp",  width = PLOT.SIZE, height = PLOT.SIZE)
                                           )
                                    )
                                  )
)