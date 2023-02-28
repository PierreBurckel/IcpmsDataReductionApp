#' Module responsible for defining how blank subtraction should be done
#' 
#' @return A list containing a reactiveVal @blankModifiers
#' 
#' blankModifiers$ @blankModifiers() (list, DataModifier)

blankProcessing_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Blank verification/processing",
          sidebarLayout(
            sidebarPanel(
              selectInput(ns("sliderInput_BlankTab_processOrView"), "Select a mode:",
                          c("View" = "view", "Process" = "process")),
              selectInput(ns("sliderInput_BlankTab_rowsToReplace"), "Replace what:",
                          "", selected = ""),
              selectInput(ns("sliderInput_BlankTab_replacementMethod"), "Replacement method:",
                          c("None" = "none", "Average blanks" = "average",
                            "Previous blank" = "previous", "Average in index" = "averageInBlankIndex")),
              selectInput(ns("sliderInput_BlankTab_rowsToReplaceFrom"), "Replace in:",
                          "All", selected = "All"),
              actionButton(ns("actionButton_BlankTab_replace"), "Set blank interpolation")),
            mainPanel(
              DT::DTOutput(ns("blankTab_table"))
            )
          )
  )
}

blankProcessing_server <- function(id, fileUpload, reactiveExpressions, indexCreation) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      blankModifiers <- reactiveVal()
      
      process <- reactive({reactiveExpressions$process})
      parameters <- reactive({fileUpload$parameters})
      customIndex <- reactive({indexCreation$custom})
      
      activeBlankModifier <- reactive({
        list(DataModifier$new(modificationMethod = input$sliderInput_BlankTab_replacementMethod,
                              modificationArguments = list(
                                linesToBeReplaced = customIndex()[[input$sliderInput_BlankTab_rowsToReplace]],
                                columnsToBeReplaced = 1:parameters()$analyteNumber,
                                linesUsedForReplacement = customIndex()[[input$sliderInput_BlankTab_rowsToReplaceFrom]],
                                columnsUsedForReplacement = 1:parameters()$analyteNumber
                                )
                              )
        )
      })
      # liveReplaceBlkTable <- reactive({
      #   
      #   process$analyteToIstdRatio()$applyModifications(activeBlankModifier())
      #   
      #   # replaceIndexWhat = customIndex()[[input$sliderInput_BlankTab_rowsToReplace]]
      #   # replaceIndexIn = customIndex()[[input$sliderInput_BlankTab_rowsToReplaceFrom]]
      #   # 
      #   # replaceValues(process$analyteToIstdRatio(), input$sliderInput_BlankTab_replacementMethod, replaceIndexWhat, replaceIndexIn, parameters())
      # })
      
      #Render ISTD table if all conditions are met
      output$blankTab_table <- DT::renderDT(datatable({
        
        blank_processOrView = input$sliderInput_BlankTab_processOrView
        
        if (blank_processOrView == "view") 
        {
          modifiedAnalyteToIstdRatio <- process()$analyteToIstdBlankRatio()
          blankTab_table <- cbind(parameters()[["categoricalDataAndTime"]][ , "Sample Name"], modifiedAnalyteToIstdRatio$getEstimation())
          names(blankTab_table) <- c("Sample Name", names(blankTab_table)[2:length(blankTab_table)])
          blankTab_table <- format(blankTab_table, digits = 3, scientific=T)
        }
        if (blank_processOrView == "process") 
        {
          activelyModifiedAnalyteToIstdRatio <- process()$analyteToIstdRatio() %>% applyModifierToEudc(activeBlankModifier())
          blankTab_table <- cbind(parameters()[["categoricalDataAndTime"]][ , "Sample Name"], activelyModifiedAnalyteToIstdRatio$getEstimation())
          names(blankTab_table) <- c("Sample Name", names(blankTab_table)[2:length(blankTab_table)])
          blankTab_table <- format(blankTab_table, digits = 3, scientific=T)
        }
        
        blankTab_table
        
      }, extensions = c('Scroller', 'Buttons'),
      options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                     scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                     buttons = c('copy', 'csv'),
                     columnDefs = list(list(width = '120px', targets = "_all")))))
      
      #Defines a proxy for changing selections in the table
      blankTab_tableProxy <- DT::dataTableProxy("blankTab_table", session = session)
      
      observe({
        input$sliderInput_BlankTab_processOrView
        input$sliderInput_BlankTab_rowsToReplace
        input$sliderInput_BlankTab_replacementMethod
        input$sliderInput_BlankTab_rowsToReplaceFrom
        if (!is.null(input$sliderInput_BlankTab_rowsToReplace)){
          DT::selectRows(blankTab_tableProxy, customIndex()[[input$sliderInput_BlankTab_rowsToReplace]])
        }
      })
      
      #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
      observeEvent(input$actionButton_BlankTab_replace, {
        blankModifiers(c(blankModifiers(), activeBlankModifier()))
        # rowReplacementIndex = customIndex()[[input$sliderInput_BlankTab_rowsToReplace]]
        # 
        # req(!is.null(liveReplaceBlkTable()) && !is.null(rowReplacementIndex))
        # 
        # blankRatioSignal <- process$analyteToIstdBlankRatio$getEstimation()
        # blankRatioRsd <- process$analyteToIstdBlankRatio$getRsd()
        # blankRatioSignal[rowReplacementIndex, ] <- liveReplaceBlkTable()$getEstimation()[rowReplacementIndex, , drop = FALSE]
        # blankRatioRsd[rowReplacementIndex, ] <- liveReplaceBlkTable()$getRsd()[rowReplacementIndex, , drop = FALSE]
        # 
        # process$analyteToIstdBlankRatio <- EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
        #                                                                        estimatedData = blankRatioSignal,
        #                                                                        uncertaintyData = blankRatioRsd,
        #                                                                        uncertaintyType = "rsd")
      })
      
      observe({
        if(input$sliderInput_BlankTab_processOrView == "process") {
          shinyjs::enable("actionButton_BlankTab_replace")
          shinyjs::enable("sliderInput_BlankTab_replacementMethod")
          shinyjs::enable("sliderInput_BlankTab_rowsToReplaceFrom")
        }
        else if(input$sliderInput_BlankTab_processOrView == "view") {
          shinyjs::disable("actionButton_BlankTab_replace")
          shinyjs::disable("sliderInput_BlankTab_replacementMethod")
          shinyjs::disable("sliderInput_BlankTab_rowsToReplaceFrom")
        }
      })
      
      observeEvent(customIndex(), {
        if (is.null(customIndex())){return()}
        updateSelectInput(session,"sliderInput_BlankTab_rowsToReplace", label  = "Replace what:", choices=names(customIndex()),names(customIndex())[1])
        updateSelectInput(session,"sliderInput_BlankTab_rowsToReplaceFrom", label  = "Replace in:", choices=names(customIndex()),names(customIndex())[1])
      })
      
      return(list(blankModifiers = blankModifiers))
    }
  )
}

