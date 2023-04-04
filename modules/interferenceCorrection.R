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
  
              actionButton(ns("actionButton_InterferenceTab_applyCorrection"), "Apply correction"),
              actionButton(ns("delete_modifier"), "Remove correction")),
            mainPanel(
              DT::DTOutput(ns("interference_modifiers_info_table"))
              )
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
      
      # Output of the information table on created index
      output$interference_modifiers_info_table <- DT::renderDT(datatable({
        interference_modifiers_info_table <- NULL
        for (modifier in interferenceModifiers()) {
          # Get the modifier information to be displayed in the table
          columnsUsedForReplacement <- modifier$getColumnsUsedForReplacement()
          columnsToBeReplaced <- modifier$getColumnsToBeReplaced()
          linesUsedForReplacement <- modifier$getLinesUsedForReplacement()
          # Convert the numerical format (line and column number) to a character format (sample and analyte names)
          columnsUsedForReplacement <- parameters()$analyteNames[columnsUsedForReplacement]
          columnsToBeReplaced <- parameters()$analyteNames[columnsToBeReplaced]
          linesUsedForReplacement <- parameters()$categoricalDataAndTime[ ,"Sample Name",drop = TRUE][linesUsedForReplacement]
          # Concatenate sample and analyte strings of length > 1 
          columnsUsedForReplacement <- paste(columnsUsedForReplacement, collapse = ", ")
          columnsToBeReplaced <- paste(columnsToBeReplaced, collapse = ", ")
          linesUsedForReplacement <- paste(linesUsedForReplacement, collapse = ", ")
          # Place the information in the table
          interference_modifiers_info_table <- rbind(interference_modifiers_info_table,
                                                     c(columnsUsedForReplacement, columnsToBeReplaced, linesUsedForReplacement))
        }
        # Rename table column names if table exists
        if (!is.null(interference_modifiers_info_table)) {
          colnames(interference_modifiers_info_table) <- c("Interfering Element", "Interfered Element", "Samples used for interference correction")
        }
        # Return the table after iteration complete
        interference_modifiers_info_table
      }, extensions = c('Scroller', 'Buttons'),
      options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                     scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                     buttons = c('copy', 'csv'),
                     columnDefs = list(list(width = '120px', targets = "_all")))))
      
      interferenceModifierTableProxy <- DT::dataTableProxy("interference_modifiers_info_table", session = session)
      
      # Handle the deletion of modifiers
      observeEvent(input$delete_modifier, {
        interferenceModifiers(interferenceModifiers()[-input$interference_modifiers_info_table_rows_selected])
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