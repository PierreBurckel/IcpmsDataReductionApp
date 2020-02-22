library(plotly)
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

ICPMS_ui <- shinyUI({
  fluidPage(
    navbarPage("ICP-MS processing",
               
               tabPanel(
                 "File upload and parameters",
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("file", "File upload"),
                     actionButton("setAsMain", "Set as main file"),
                     actionButton("setAsStd", "Set as standard file"),
                     downloadButton("downloadISTDTemplate", "Download ISTD template"),
                     actionButton("setAsISTD", "Set as ISTD file"),
                     actionButton("extract", "Extract"),
                     htmlOutput("extract_ready_txt"),
                     htmlOutput("ISTD_not_extracted_txt")),
                   mainPanel(
                     htmlOutput("main_assignment_txt"),
                     htmlOutput("std_assignment_txt"),
                     htmlOutput("ISTD_assignment_txt"),
                     div(style="display:inline-block",numericInput(inputId="fileUpload_nrow", label="Rows", value = 10, min = 1)),
                     div(style="display:inline-block",numericInput(inputId="fileUpload_ncolumn", label="Columns", value = 6, min = 1)),
                     tableOutput("tempFilePreview")
                   )
                 )
               ),
               
               tabPanel("Index creation",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("searchIndexwhere", "In:", c("Sample Names" = "smp", "Levels" = "lvl", "Type" = "type")),
                            textInput("searchIndexwhat", "Search:", ""),
                            selectInput("searchIndexhow", "Using:", c("Regular expression" = "regexp", "Exact match" = "ematch")),
                            textInput("searchIndexName", "Index Name:", ""),
                            actionButton("indexSelectAll", "Select All"),
                            actionButton("searchIndexCreate", "Create new index")),
                          mainPanel(
                            selectInput("searchIndexDisplay", "Show:", c("Internal standards" = "ISTD", "Analytes" = "analytes")),
                            DT::DTOutput("index_table")
                          )
                        )
               ),
               
               tabPanel("ISTD verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            shinyjs::useShinyjs(),
                            sliderInput("ISTDcolSlider", label = NULL, min = 1, 
                                      max = 100, value = c(1, 1)),
                            selectInput("ISTDinteractionMode", "Select a mode:",
                                        c("View" = "view", "Process" = "process")),
                            selectInput("indexISTDchoiceWhat", "Replace what:",
                                        "", selected = ""),
                            selectInput("ISTDcorrectionMethod", "Replace calibration ISTD:",
                                        c("None" = "none", "Average blanks" = "mean",
                                          "Previous blank" = "previous")),
                            selectInput("indexISTDchoiceIn", "Replace in:",
                                        "All", selected = "All"),
                            actionButton("setISTDcorrectionMethod", "Set ISTD correction"),
                            checkboxGroupInput("ISTDmodifiers", "Modifiers:", choices = NULL, selected = NULL)),
                          mainPanel(
                            DT::DTOutput("ISTDtable")
                          )
                        )
               ),
               
               tabPanel("Blank verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            shinyjs::useShinyjs(),
                            checkboxInput("useBlankCorrection", "Check to perform blank correction", value=FALSE),
                            selectInput("blkInteractionMode", "Select a mode:",
                                        c("View" = "view", "Process" = "process")),
                            selectInput("indexBlkchoiceWhat", "Replace what:",
                                        "", selected = ""),
                            selectInput("blkInterpolationMethod", "Interpolate blanks:",
                                        c("None" = "none", "Average blanks" = "mean",
                                          "Previous blank" = "previous")),
                            selectInput("indexBlkchoiceIn", "Replace in:",
                                        "All", selected = "All"),
                            actionButton("setBlkInterpolationMethod", "Set blank interpolation"),
                            checkboxGroupInput("blankModifiers", "Modifiers:", choices = NULL, selected = NULL)),
                          mainPanel(
                            div(style="display:inline-block",sliderInput("blkColSlider", label = h3("Column range"), min = 1, 
                                                                         max = 100, value = c(1, 1))),
                            DT::DTOutput("blkTable")
                          )
                        )
               ),
               
               tabPanel("Calibration verification",
                        sidebarLayout(
                          sidebarPanel(
                            checkboxInput("forceIntercept", "Force intercept", value = FALSE),
                            checkboxInput("autoAdaptCalibration", "Use min/max standards", value = FALSE),
                            checkboxInput("useWeithedRegression", "Use weighted linear regression", value = FALSE),
                            selectInput("regressionWeight", "Weight:", c("1/var" = "1/var", "1/Y" = "1/Y", "1/max(var,Y)" = "1/max(var,Y)"), selected = "1/var"),
                            actionButton("setRegression", "Set regression"),
                            actionButton("setRegressionALL", "Set to all")),
                          mainPanel(
                            selectInput("calibrationElement", "Element:", "", selected = ""),
                            plotlyOutput("calibrationPlot")
                          )
                        )
               ),
               
               tabPanel("Drift verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("selectDriftIndex", "Define drift index:", "", selected = ""),
                            actionButton("setAsDriftIndex", "Set drift index"),
                            selectInput("e_drift", "Choose element:", "", selected = ""),
                            numericInput("e_ind_drift", "Element number", 1),
                            selectInput("driftModelSelection", "Choose regression:", c("None" = "None", "Linear" = "Linear", "Quadratic" = "Quadratic"), selected = "None"),
                            actionButton("setDriftCorrection", "Set ISTD correction"),
                            textOutput("warningDrifr")),
                          mainPanel(
                            tableOutput("smpBlkCorTable"),
                            plotlyOutput("driftPlot")
                          )
                        )
               ),
               
               tabPanel("Process",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("viewConcentrationSwitch", "Display:",
                                        c("Concentrations (ppb)" = 1, "RSD (%)" = 2, "Columnwise concatenation" = 3)),
                            selectInput("viewConcentrationIndex", "View index:",
                                        "All", selected = "All"),
                            actionButton("process", "Process data"),
                            downloadButton("downloadData", "Download table")),
                          mainPanel(
                            tableOutput("conc")
                          )
                        )
               )
    )
  )
})