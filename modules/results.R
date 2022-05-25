results_ui <- function(id) {
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
}

results_server <- function(input, output, session) {
  
}