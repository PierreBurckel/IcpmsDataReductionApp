#' Module responsible for formatting, viewing and downloading final concentration
#' and intermediate results

results_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Process",
          sidebarLayout(
            sidebarPanel(
              selectInput(ns("viewCustomDataSwitch"), "Display:",
                          c("Value" = 1, "Uncertainty" = 2, "Columnwise concatenation" = 3)),
              selectInput(ns("viewConcentrationIndex"), "View index:",
                          "All", selected = "All"),
              selectInput(ns("dataReductionStateSelection"), "Select the type of data to download:",
                          c("Counts per second",
                            "Interference counts per second", "Interference corrected counts per second",
                            "Internal Standard matrix (adapted for analytes)",
                            "Internal Standard ratio", "Blank matrix", "Blank corrected signal (CPS or ratio)",
                            "Detection Limit", "Quantification Limit",
                            "Drift matrix", "Drift corrected signal (CPS or ratio, blank corrected)", "Concentration"),
                          selected = "Concentration"),
              downloadButton(ns("downloadCustomTable"), "Download custom table")),
            mainPanel(
              DT::DTOutput(ns("customDataTable"))
            )
          )
  )
}

results_server <- function(id, indexCreation, reactiveExpressions, fileUpload) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      customIndex <- reactive({indexCreation$custom})
      process <- reactive({reactiveExpressions$process})
      parameters <- reactive({fileUpload$parameters})
      
      downloadTabData <- reactive({
        
        selectedIndex <- customIndex()[[input$viewConcentrationIndex]]
        
        downloadTabDataSignal <- switch(input$dataReductionStateSelection,
                                        "Counts per second" = process()$analyteCountsPerSecond(),
                                        "Interference counts per second" = process()$interferenceCountsPerSecond()$getEstimation(),
                                        "Interference corrected counts per second" = process()$interferenceCorrectedCountsPerSecond()$getEstimation(),
                                        "Internal Standard matrix (adapted for analytes)" = process()$internalStandardEudcAdaptedToAnalytes()$getEstimation(),
                                        "Internal Standard ratio" = process()$analyteToIstdRatio()$getEstimation(),
                                        "Blank matrix" = process()$analyteToIstdBlankRatio()$getEstimation(),
                                        "Detection Limit" = process()$detectionLimit()$getEstimation(),
                                        "Quantification Limit" = process()$quantificationLimit()$getEstimation(),
                                        "Blank corrected signal (CPS or ratio)" = process()$analyteToIstdRatioBlankCorrected()$getEstimation(),
                                        "Drift matrix" = process()$driftFactorEudc()$getEstimation(),
                                        "Drift corrected signal (CPS or ratio, blank corrected)" = process()$driftAndBlankCorrectedEudc()$getEstimation(),
                                        "Concentration" = process()$concentration()$getEstimation()
        )
        
        downloadTabDataSignal[downloadTabDataSignal < 0] <- "<blk"
        downloadTabDataSignal[is.na(downloadTabDataSignal)] <- "N/A"
        
        downloadTabDataRsd <- switch(input$dataReductionStateSelection,
                                     "Counts per second" = process()$analyteCountsPerSecondRelativeStandardDeviation(),
                                     "Interference counts per second" = process()$interferenceCountsPerSecond()$getRsd(),
                                     "Interference corrected counts per second" = process()$interferenceCorrectedCountsPerSecond()$getRsd(),
                                     "Internal Standard matrix (adapted for analytes)" = process()$internalStandardEudcAdaptedToAnalytes()$getRsd(),
                                     "Internal Standard ratio" = process()$analyteToIstdRatio()$getRsd(),
                                     "Blank matrix" = process()$analyteToIstdBlankRatio()$getRsd(),
                                     "Detection Limit" = process()$detectionLimit()$getRsd(),
                                     "Quantification Limit" = process()$quantificationLimit()$getRsd(),
                                     "Blank corrected signal (CPS or ratio)" = process()$analyteToIstdRatioBlankCorrected()$getRsd(),
                                     "Drift matrix" = process()$driftFactorEudc()$getRsd(),
                                     "Drift corrected signal (CPS or ratio, blank corrected)" = process()$driftAndBlankCorrectedEudc()$getRsd(),
                                     "Concentration" = process()$concentration()$getRsd()
        )
        
        downloadTabDataRsd[downloadTabDataSignal < 0] <- "N/A"
        downloadTabDataRsd[is.na(downloadTabDataRsd)] <- "N/A"
        
        colnames(downloadTabDataSignal) <- parameters()[["analyteNames"]]
        colnames(downloadTabDataRsd) <- parameters()[["analyteNames"]]
        
        mergedSignalAndRsd <- mergeMatrixes(downloadTabDataSignal, downloadTabDataRsd)
        mergedSignalAndRsdColumnNames <- mergeMatrixes(matrix(parameters()[["analyteNames"]], nrow = 1, ncol = parameters()[["analyteNumber"]]),
                                                       matrix(paste0(parameters()[["analyteNames"]], " RSD (%)"), nrow = 1, ncol = parameters()[["analyteNumber"]]))
        colnames(mergedSignalAndRsd) <- mergedSignalAndRsdColumnNames
        
        
        customData <- switch(input$viewCustomDataSwitch,
                             "1" = downloadTabDataSignal,
                             "2" = downloadTabDataRsd,
                             "3" = mergedSignalAndRsd
        )
        
        customDataWithHeader <- cbind(parameters()[["categoricalDataAndTime"]][, "Sample Name", drop = FALSE], customData)
        
        return(customDataWithHeader[selectedIndex, ])
        
      })
      
      output$customDataTable <- DT::renderDT(datatable({
        
        downloadTabData()
        
      }, extensions = c('Scroller', 'Buttons'),
      options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                     scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                     buttons = c('copy', 'csv'),
                     columnDefs = list(list(width = '120px', targets = "_all")))))
      
      
      # Downloadable data -------------------------------------------------------
      
      output$downloadCustomTable <- downloadHandler(
        
        filename =  paste('data_', input$viewConcentrationIndex, '_', Sys.Date(), '.csv', sep=''),
        content = function(file) 
        {
          write.table(downloadTabData(), file, sep=parameters()$csvDelimitation, quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
      )
      
      observeEvent(customIndex(), {
        if (is.null(customIndex())){return()}
        updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(customIndex())),"All")
      })
    }
  )
}