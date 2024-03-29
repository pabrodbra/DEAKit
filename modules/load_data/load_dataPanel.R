load.data.panel <- 
  tabPanel(title = "Load Data",
           fluidRow(
             column(12,
                    span("Load Data parameters:"),
                     wellPanel(
                       div(class = "col-md-12",
                       #div(class = "col-md-6",
                           
                       
                         # 1 - RCC Data directory
                         textInput("rcc.path", "RCC Directory:", value = file.path(getwd(), "data/PanCancer-Pathways-57"),
                                   placeholder = "/home/userXY/data/RCC_files/"),
                         helpText("Path of the Directory for the RCC files. Names must be identified as SAMPLE[_keyX][_keyX+1][...].RCC"),
                         
                         # 2 - RLF Filename
                         textInput("rlf.path", "RLF File:", value = file.path(getwd(), "data/PanCancer-Pathways-57/NS_CancerPath_C2535.rlf"),
                                   placeholder = "/home/userXY/data/RLF_file.rlf"),
                         helpText("Path of the Location of the RLF File"),
                         
                         # 3 - KEY OF INTEREST - For all
                         numericInput("key.of.interest", "Key", value = 1, min=1),
                         helpText("Key number to compare")
                         ,
                       #),
                       
                       #div(class = "col-md-6", 
                         # 4 - LABEL OF SUCH KEY - For all
                         textInput("label.of.key", "Label for Key", value = "PCR",
                                   placeholder = "PCR"),
                         helpText("Label for such key"),
                         
                         # 5 - VALUES OF SUCH KEY TO COMPARE - For all 
                         textInput("key.value.1", "Possible condition 1 of Key X", value = 0,
                                   placeholder = "0"),
                         helpText("Condition 1 for Selected Key"),
                         
                         textInput("key.value.2", "Possible condition 2 of Key X", value = 1,
                                   placeholder = "1"),
                         helpText("Condition 2 for Selected Key")
                         
                       ),
                       
                       div(class="text-center", 
                           useShinyjs(),
                           tags$style(appCSS),
                           withBusyIndicatorUI(
                             actionButton("loadDataButton", 
                                          "Load Data", 
                                          class="btn-primary"))
                       )
                     )
             ),
             column(12,
                    span("Load Data results:"),
                    wellPanel(id = "ld.res",
                      span("Metadata:"),
                      dataTableOutput("loaded.metadata")),
                    wellPanel(
                      span("Count matrix (first 5 columns) :"),
                      dataTableOutput("loaded.counts")
                    )
             )
           )
  )
