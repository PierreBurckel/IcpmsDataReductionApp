interferenceCorrection_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Interference correction",
          sidebarLayout(
            sidebarPanel(
              selectInput(ns("sliderInput_InterferenceTab_indexForCorrectionFactorCalculation"), "Compute correction factor from:",
                          "", selected = ""),
              selectInput(ns("sliderInput_InterferenceTab_interferingElement"), "Select interfering element:",
                          "", selected = ""),
              selectInput(ns("sliderInput_InterferenceTab_interferedElement"), "Select interfered element:",
                          "", selected = ""),
  
              actionButton(ns("actionButton_InterferenceTab_applyCorrection"), "Apply correction")),
            mainPanel()
          )
  )
}

interferenceCorrection_server <- function(id, fileUpload, reactiveExpressions, indexCreation) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      interferenceModifiers <- reactiveVal()
      
      process <- reactive({reactiveExpressions$process})
      parameters <- reactive({fileUpload$parameters})
      customIndex <- reactive({indexCreation$custom})
      
      activeInterferenceModifier <- reactive({
        list(
          DataModifier$new(modificationMethod = "interferenceCorrection",
                              modificationArguments = list(
                                linesToBeReplaced = 1:parameters()$sampleNumber,
                                columnsToBeReplaced = which(input$sliderInput_InterferenceTab_interferedElement == parameters()$analyteNames),
                                linesUsedForReplacement = customIndex()[[input$sliderInput_InterferenceTab_indexForCorrectionFactorCalculation]],
                                columnsUsedForReplacement = which(input$sliderInput_InterferenceTab_interferingElement == parameters()$analyteNames)
                              )
          )
        )
      })
      
      observeEvent(input$actionButton_InterferenceTab_applyCorrection, {
        # Merging the current DataModifier list in the interferenceModifiers list
        interferenceModifiers(c(interferenceModifiers(), activeInterferenceModifier()))
      })
      
      observeEvent(customIndex(), {
        if (is.null(customIndex())){return()}
        updateSelectInput(session,"sliderInput_InterferenceTab_indexForCorrectionFactorCalculation",
                          label  = "Compute correction factor from:", choices=names(customIndex()),names(customIndex())[1])
      })
      
      observeEvent(parameters()$analyteNames, {
        if (is.null(parameters())){return()}
        updateSelectInput(session,"sliderInput_InterferenceTab_interferingElement",
                          label  = "Select interfering element:",
                          choices=parameters()$analyteNames,
                          parameters()$analyteNames[1])
        updateSelectInput(session,"sliderInput_InterferenceTab_interferedElement",
                          label  = "Select interfered element:",
                          choices=parameters()$analyteNames,
                          parameters()$analyteNames[1])
      })
      
      return(list(interferenceModifiers = interferenceModifiers))
    }
  )
  
}