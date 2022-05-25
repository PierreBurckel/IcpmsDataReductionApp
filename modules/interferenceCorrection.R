interferenceCorrection_ui <- function(id) {
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
  )
}

interferenceCorrection_server <- function(input, output, session) {
  
}