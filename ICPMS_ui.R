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
                     htmlOutput("extractionReadyText")),
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
               
               tabPanel("Blank verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("blkInteractionMode", "Select a mode:",
                                        c("View" = "view", "Process" = "process")),
                            selectInput("indexBlkchoiceWhat", "Replace what:",
                                        "", selected = ""),
                            selectInput("blkInterpolationMethod", "Replacement method:",
                                        c("None" = "none", "Average blanks" = "mean",
                                          "Previous blank" = "prev", "Average in index" = "averageInBlankIndex")),
                            selectInput("indexBlkchoiceIn", "Replace in:",
                                        "All", selected = "All"),
                            actionButton("setBlkInterpolationMethod", "Set blank interpolation")),
                          mainPanel(
                            div(style="display:inline-block",sliderInput("blkColSlider", label = h3("Column range"), min = 1, 
                                                                         max = 100, value = c(1, 1))),
                            DT::DTOutput("blkTable")
                          )
                        )
               ),
               
               tabPanel("Drift verification/processing",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("selectDriftIndex", "Define drift index:", "", selected = ""),
                            actionButton("setAsDriftIndex", "Set drift index"),
                            selectInput("e_drift", "Choose element:", "", selected = ""),
                            actionButton("changeElementDisplayedForDrift", "Switch to element"),
                            numericInput("e_ind_drift", "Element number", 1),
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
                                        c("Counts per second", "Internal Standard matrix (adapted for analytes)",
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