driftProcessing_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Drift verification/processing",
           sidebarLayout(
             sidebarPanel(
               selectInput(ns("selectDriftIndex"), "Define drift index:", "", selected = ""),
               actionButton(ns("setAsDriftIndex"), "Set drift index"),
               selectInput(ns("driftTab_selectInput_analyteName"), "Choose element:", "", selected = ""),
               actionButton(ns("changeElementDisplayedForDrift"), "Switch to element"),
               numericInput(ns("driftTab_numericInput_analyteNumber"), "Element number", 1),
               actionButton(ns("setDriftCorrection"), "Set ISTD correction"),
               textOutput(ns("warningDrifr")),
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
               DT::DTOutput(ns("smpBlkCorTable")),
               plotOutput(ns("driftPlot")),
               tableOutput(ns("test"))
             )
           )
  )
}

driftProcessing_server <- function(id, fileUpload, indexCreation, reactiveExpressions) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      browser()
      customIndex <- reactive({indexCreation$custom})
      driftIndex <- reactive({indexCreation$drift})
      process = reactiveExpressions$process
      parameters = fileUpload$parameters
      
      observeEvent(input$setAsDriftIndex, {
        if (input$selectDriftIndex != ""){
          driftIndex() <- customIndex()[[input$selectDriftIndex]]
        }
      })
      
      output$smpBlkCorTable <- DT::renderDT(datatable({
        if (is.null(process$analyteToIstdRatioBlankCorrected()$getEstimation()) | is.null(driftIndex()) | !is.integer(input$driftTab_numericInput_analyteNumber)){return()}
        driftSignal <- process$analyteToIstdRatioBlankCorrected()$getEstimation()
        if (input$driftTab_numericInput_analyteNumber < 5){
          lc <- 1
        }
        else{
          lc <- max(which((1:input$driftTab_numericInput_analyteNumber)%%5 == 0))
        }
        uc <- min(lc + 4, parameters$analyteNumber)
        driftTable <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"][driftIndex()],driftSignal[driftIndex(),lc:uc])
        names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
        driftTable
        
      }, extensions = c('Scroller', 'Buttons'),
      options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                     scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                     buttons = c('copy', 'csv'),
                     columnDefs = list(list(width = '120px', targets = "_all")))))
      
      
      observeEvent(parameters$analyteNames, {
        if (is.null(parameters$analyteNames)){return()}
        updateSelectInput(session,"sliderInput_InterferenceTab_interferingElement",choices=parameters$analyteNames,selected=parameters$analyteNames[1])
        updateSelectInput(session,"sliderInput_InterferenceTab_interferedElement",choices=parameters$analyteNames,selected=parameters$analyteNames[1])
        updateSelectInput(session,"driftTab_selectInput_analyteName",choices=parameters$analyteNames,selected=parameters$analyteNames[1])
        updateNumericInput(session,"driftTab_numericInput_analyteNumber", "Element number", 1, min = 1, max = length(parameters$analyteNames))
      })
      
      observeEvent(input$changeElementDisplayedForDrift, {
        if (is.null(input$driftTab_selectInput_analyteName)){return()}
        updateNumericInput(session, "driftTab_numericInput_analyteNumber", value = grep(input$driftTab_selectInput_analyteName,parameters$analyteNames,fixed=TRUE))
      })
      
      observeEvent(input$driftTab_numericInput_analyteNumber, {
        if (is.null(input$driftTab_numericInput_analyteNumber)){return()}
        updateSelectInput(session,"driftTab_selectInput_analyteName",selected=parameters$analyteNames[input$driftTab_numericInput_analyteNumber])
      })
      
      yplus <- reactive({(process$analyteToIstdRatioBlankCorrected()$getEstimation() + process$analyteToIstdRatioBlankCorrected()$getSd())[driftIndex(), input$driftTab_numericInput_analyteNumber]})
      yminus <- reactive({(process$analyteToIstdRatioBlankCorrected()$getEstimation() - process$analyteToIstdRatioBlankCorrected()$getSd())[driftIndex(), input$driftTab_numericInput_analyteNumber]})
      
      output$driftPlot <- renderPlot({
        if (is.null(process$analyteToIstdRatioBlankCorrected()$getEstimation()) | is.null(driftIndex())){return()}
        
        elementFullName <- parameters$analyteNames[input$driftTab_numericInput_analyteNumber]
        
        driftTime <- parameters$deltaTime()[driftIndex()]
        driftValue <- process$analyteToIstdRatioBlankCorrected()$getEstimation()[driftIndex(),input$driftTab_numericInput_analyteNumber]
        driftSd <- process$analyteToIstdRatioBlankCorrected()$getSd()[driftIndex(),input$driftTab_numericInput_analyteNumber]
        
        elementSpecificDriftIndex <- parameters$listOfElementSpecificDriftIndex()[[elementFullName]]
        timeOfDriftAfterFirstStandard <- parameters$deltaTime()[elementSpecificDriftIndex]
        valueOfDriftAfterFirstStandard <- process$analyteToIstdRatioBlankCorrected()$getEstimation()[elementSpecificDriftIndex, input$driftTab_numericInput_analyteNumber]
        
        colors <- c("Drift values" = "lightblue", "Drift values after first standard" = "red")
        
        driftPlot <- ggplot() + xlab("Time") + ylab(elementFullName)
        driftPlot <- driftPlot + 
          geom_errorbar(data = data.frame(driftTime, driftValue, driftSd),
                        aes(driftTime, driftValue, ymin=driftValue - driftSd, ymax=driftValue + driftSd), width=.1) + 
          geom_point(data = data.frame(driftTime, driftValue),
                     aes(driftTime, driftValue, color = "Drift values"), size = 10) + 
          
          geom_point(data = data.frame(timeOfDriftAfterFirstStandard, valueOfDriftAfterFirstStandard),
                     aes(timeOfDriftAfterFirstStandard, valueOfDriftAfterFirstStandard, color = "Drift values after first standard"), size = 6)
        
        if (parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber] == TRUE && class(process$elementSpecificDriftModels()[[elementFullName]]) == "lm"){
          time_0 = driftTime[1]
          time_f = driftTime[length(driftTime)]
          timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
          
          dt = as.numeric(driftTime)
          driftModel <- process$elementSpecificDriftModels()[[elementFullName]]
          
          driftPredict=predict(driftModel, newdata = data.frame(dt = timeInterval))
          
          driftPlot <- driftPlot + geom_line(data = data.frame(timeInterval, driftPredict),
                                             aes(timeInterval, driftPredict), color = "red")
        }
        driftPlot + scale_color_manual(values = colors)
      })
      
      observeEvent(input$setDriftCorrection, {
        if (is.null(parameters$driftCorrectedElements)){return()}
        parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber] <- !parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber]
      })
      
      observeEvent(input$pressedKey, {
        if (is.null(parameters$driftCorrectedElements) || input$pressedKeyId != C_LETTER_KEYCODE || input$tagId != "driftTab_numericInput_analyteNumber"){return()}
        parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber] <- !parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber]
      })
      
    }
  )
  
  
}

