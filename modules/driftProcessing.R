driftProcessing_ui <- function(id) {
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
  )
}

driftProcessing_server <- function(input, output,server) {
  
}

