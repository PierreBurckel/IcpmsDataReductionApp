# library(plotly)
# source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

ICPMS_ui <- shinyUI({
  fluidPage(
    useShinyalert(),
    shinyjs::useShinyjs(),
    navbarPage("ICP-MS processing",
               
               tabPanel(
                 "File upload and parameters",
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("uploadedFile", "File upload"),
                     actionButton("setAsMainFile", "Set as main file"),
                     actionButton("setAsStandardFile", "Set as standard file"),
                     downloadButton("downloadISTDTemplate", "Download ISTD template"),
                     actionButton("setAsInternalStandardFile", "Set as ISTD file"),
                     actionButton("extract", "Extract"),
                     htmlOutput("extractionReadyText"),
                     selectInput(inputId = "csvDelimitation", label = "Select csv field separator (the separator will be used during extraction,
                                 you therefore need to use the same separator for all your files)",
                                 choices = c(";", ","), selected = ";"),
                     textInput(inputId = "dateFormat", label = "Date format", value = "%d/%m/%Y %H:%M")),
                   mainPanel(
                     htmlOutput("mainFileAssignmentText"),
                     htmlOutput("standardFileAssignmentText"),
                     htmlOutput("internalStandardFileAssignmentText"),
                     div(style="display:inline-block", numericInput(inputId="fileUpload_nrow", label="Rows", value = 10, min = 1)),
                     div(style="display:inline-block", numericInput(inputId="fileUpload_ncolumn", label="Columns", value = 6, min = 1)),
                     tableOutput("uploadedFilePreviewTable")
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
                            DT::DTOutput("indexTable")
                          )
                        )
               ),
               
               tabPanel("Interference correction",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("sliderInput_InterferenceTab_indexForCorrectionFactorCalculation", "Compute correction factor from:",
                                        "", selected = ""),
                            selectInput("sliderInput_InterferenceTab_interferingElement", "Select interfering element:",
                                        "", selected = ""),
                            selectInput("sliderInput_InterferenceTab_interferedElement", "Select interfered element:",
                                        "", selected = ""),
          
                            actionButton("actionButton_InterferenceTab_applyCorrection", "Apply correction")),
                          mainPanel()
                        )
               ),
               
               tabPanel("Blank verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("sliderInput_BlankTab_processOrView", "Select a mode:",
                                        c("View" = "view", "Process" = "process")),
                            selectInput("sliderInput_BlankTab_rowsToReplace", "Replace what:",
                                        "", selected = ""),
                            selectInput("sliderInput_BlankTab_replacementMethod", "Replacement method:",
                                        c("None" = "none", "Average blanks" = "mean",
                                          "Previous blank" = "prev", "Average in index" = "averageInBlankIndex")),
                            selectInput("sliderInput_BlankTab_rowsToReplaceFrom", "Replace in:",
                                        "All", selected = "All"),
                            actionButton("actionButton_BlankTab_replace", "Set blank interpolation")),
                          mainPanel(
                            checkboxInput("enableBlankCorrection", "Enable blank correction", value = TRUE),
                            DT::DTOutput("blankTab_table")
                          )
                        )
               ),
               
               tabPanel("Drift verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("selectDriftIndex", "Define drift index:", "", selected = ""),
                            actionButton("setAsDriftIndex", "Set drift index"),
                            selectInput("driftTab_selectInput_analyteName", "Choose element:", "", selected = ""),
                            actionButton("changeElementDisplayedForDrift", "Switch to element"),
                            numericInput("driftTab_numericInput_analyteNumber", "Element number", 1),
                            actionButton("setDriftCorrection", "Set ISTD correction"),
                            textOutput("warningDrifr"),
                            tags$script('
                              pressedKeyCount = 0;
                              $(document).on("keydown", function (e) {
                                 var tag = e.target.id;
                                 Shiny.onInputChange("tagId", tag);
                                 Shiny.onInputChange("pressedKey", pressedKeyCount++);
                                 Shiny.onInputChange("pressedKeyId", e.which);
                              });'
                            )),
                          mainPanel(
                            DT::DTOutput("smpBlkCorTable"),
                            plotOutput("driftPlot"),
                            tableOutput("test")
                          )
                        )
               ),
               
               tabPanel("Process",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("viewCustomDataSwitch", "Display:",
                                        c("Value" = 1, "Uncertainty" = 2, "Columnwise concatenation" = 3)),
                            selectInput("viewConcentrationIndex", "View index:",
                                        "All", selected = "All"),
                            selectInput("dataReductionStateSelection", "Select the type of data to download:",
                                        c("Counts per second",
                                          "Interference counts per second", "Interference corrected counts per second",
                                          "Internal Standard matrix (adapted for analytes)",
                                          "Internal Standard ratio", "Blank matrix", "Blank corrected signal (CPS or ratio)",
                                          "Drift matrix", "Drift corrected signal (CPS or ratio, blank corrected)", "Concentration"),
                                        selected = "Concentration"),
                            downloadButton("downloadCustomTable", "Download custom table")),
                          mainPanel(
                            DT::DTOutput("customDataTable")
                          )
                        )
               )
    )
  )
})