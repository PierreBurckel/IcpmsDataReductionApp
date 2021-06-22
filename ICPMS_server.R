
C_LETTER_KEYCODE <- 67

ICPMS_server <- function(input, output, session) {
  
  shinyjs::disable(selector = '.navbar-nav a[data-value="Index creation"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Process"')
  
  #Setting up lists of reactiveValues variables
  #uploadedFile contains the information about user-uploaded files
  uploadedFile <- reactiveValues()
  #extracted$data contains the list of important data after extraction
  extracted <- reactiveValues()
  #index contains all indexes serving as masks to display/process specific lines or columns of the raw data
  rowIndexInMain <- reactiveValues()
  columnIndexInMain <- reactiveValues()
  rawData <- reactiveValues()
  process <- reactiveValues()
  parameters <- reactiveValues()
  applicationState <- reactiveValues(isExtractionSuccessful = FALSE)
  #isExtractionReady is TRUE or FALSE, depending on whether extraction can be done or not (required files assigned, variables names inputed)
  
  isValidISTD <- reactive({(!is.null(columnIndexInMain$internalStandardCountsPerSecond)) & (!is.null(extracted$internalStandard))})
  
  # columnIndexInMain$analytesCountsPerSecond <- extracted$secondRowOfMain == "CPS" & !grepl("ISTD", extracted$firstRowOfMain)
  # columnIndexInMain$analytesCountsPerSecondRelativeStandardDeviation <- extracted$secondRowOfMain == "CPS RSD" & !grepl("ISTD", extracted$firstRowOfMain)
  # columnIndexInMain$internalStandardCountsPerSecond <- extracted$secondRowOfMain == "CPS" & grepl("ISTD", extracted$firstRowOfMain)
  # columnIndexInMain$internalStandardCountsPerSecondRelativeStandardDeviation <- extracted$secondRowOfMain == "CPS RSD" & grepl("ISTD", extracted$firstRowOfMain)
  # columnIndexInMain$numericalColumns <- min(which(columnIndexInMain$analytesCountsPerSecond)):length(extracted$main)
  
  process$countsPerSecond <- reactive({
    extracted$main[ , extracted$secondRowOfMain == "CPS", drop=FALSE]
  })
  
  process$analyteCountsPerSecond <- reactive({
    process$countsPerSecond()[ , parameters$analyteNames, drop=FALSE]
  })
  
  process$internalStandardCountsPerSecond <- reactive({
    process$countsPerSecond()[ , parameters$internalStandardNames, drop=FALSE]
  })
  
  process$relativeStandardDeviation <- reactive({
    extracted$main[ , extracted$secondRowOfMain == "CPS RSD", drop=FALSE]
  })
  
  process$analyteCountsPerSecondRelativeStandardDeviation <- reactive({
    process$relativeStandardDeviation()[ , parameters$analyteNames, drop=FALSE]
  })
  
  process$internalStandardCountsPerSecondRelativeStandardDeviation <- reactive({
    process$relativeStandardDeviation()[ , parameters$internalStandardNames, drop=FALSE]
  })
  
  # elementNames <- reactive({
  #   header_1 = extracted$firstRowOfMain
  #   header_1[columnIndexInMain$analytesCountsPerSecond]})
  # 
  # elementNumber <- reactive({
  #   length(which(columnIndexInMain$analytesCountsPerSecond))})
  # 
  # ISNames <- reactive({
  #   header_1 = extracted$firstRowOfMain
  #   header_1[columnIndexInMain$internalStandardCountsPerSecond]})
  # 
  # ISNumber <- reactive({
  #   length(which(columnIndexInMain$internalStandardCountsPerSecond))})
  # 
  # sampleNumber <- reactive({
  #   nrow(extracted$main)})
  # 
  # #Important columns
  # timeColumn <- reactive({
  #   rawData = extracted$main
  #   header_2 = extracted$secondRowOfMain
  #   as.POSIXct(rawData[,which(header_2=="Acq. Date-Time")], format="%d/%m/%Y %H:%M")})
  # 
  # nameColumn <- reactive({
  #   rawData = extracted$main
  #   header_2 = extracted$secondRowOfMain
  #   rawData[,which(header_2=="Sample Name")]})
  # 
  # typeColumn <- reactive({
  #   rawData = extracted$main
  #   header_2 = extracted$secondRowOfMain
  #   rawData[,which(header_2=="Type")]})
  # 
  # levelColumn <- reactive({
  #   rawData = extracted$main
  #   header_2 = extracted$secondRowOfMain
  #   rawData[,which(header_2=="Level")]})
  # 
  # dtimeColumn <- reactive({
  #   parameters[["categoricalDataAndTime"]][ , "Time"] - parameters[["categoricalDataAndTime"]][ , "Time"][1]})
  
  internalStandardMatrixAdaptedToAnalytes <- reactive({
    if (!is.null(extracted$internalStandard)) 
    {
      return(createInternalStandardMatrixAdaptedToAnalytes(extracted$internalStandard, list(signal=process$internalStandardCountsPerSecond(), RSD=process$internalStandardCountsPerSecondRelativeStandardDeviation())))
    } 
    else
    {
      return(list(signal = matrix(1, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                  RSD = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber)))
    }
  })
  
  process$ratio <- reactive({propagateUncertainty(a = list(signal=process$analyteCountsPerSecond(), RSD=process$analyteCountsPerSecondRelativeStandardDeviation()), b = internalStandardMatrixAdaptedToAnalytes(), operation="division")})
  
  liveReplaceBlkTable <- reactive({
    
    replaceIndexWhat = rowIndexInMain$custom[[input$indexBlkchoiceWhat]]
    replaceIndexIn = rowIndexInMain$custom[[input$indexBlkchoiceIn]]
    
    replaceValues(process$ratio(), 1:parameters$analyteNumber, input$blkInterpolationMethod, replaceIndexWhat, replaceIndexIn)
    })
    
  
  process$ratio_cor_b <- reactive({
      propagateUncertainty(a=process$ratio(), b=process$blk_ratio, operation="substraction")
  })
  
  yplus <- reactive({(process$ratio_cor_b()[[1]] + process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[rowIndexInMain$drift, input$e_ind_drift]})
  yminus <- reactive({(process$ratio_cor_b()[[1]] - process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[rowIndexInMain$drift, input$e_ind_drift]})
  
  process$driftFactorMatrix <- reactive({
  
    for (elementIndex in 1:parameters$analyteNumber){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      driftIndexAfterFirstStandard <- getElementDriftIndex(elementFullName = elementFullName, stdDataFrame = extracted$standard, 
                                                           stdIdentificationColumn=parameters[["categoricalDataAndTime"]][ , "Level"], driftIndex = rowIndexInMain$drift)
      
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- parameters[["categoricalDataAndTime"]][ , "Time"] - parameters[["categoricalDataAndTime"]][ , "Time"][driftIndexAfterFirstStandard][1]
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[deltaTimeWithFirstDriftAfterFirstStandardAsOrigin < 0] <- 0
      
      driftSignalAndDeltaTime <- cbind(Signal=process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex],
                                       dt=deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[driftIndexAfterFirstStandard])  
      
      if (all(is.na(process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex]))){
        driftPredict = rep(NA, parameters$sampleNumber)
      }
      else{
        driftModel <- lm(Signal ~ poly(dt, degree=2, raw=TRUE), data = as.data.frame(driftSignalAndDeltaTime))
        driftPredict=predict(driftModel, newdata = data.frame(dt=deltaTimeWithFirstDriftAfterFirstStandardAsOrigin))
      }
      
      if (elementIndex == 1){
        driftDataFrame <- data.frame(driftPredict)
      }
      else{
        driftDataFrame[elementIndex] <- driftPredict
      }
  
      if (parameters$driftCorrectedElements[elementIndex] == TRUE){
        driftDataFrame[elementIndex] <- driftDataFrame[elementIndex] / summary(driftModel)$coefficients[1,1]
        driftDataFrame[is.na(driftDataFrame[,elementIndex]), elementIndex] <- 1
      } 
      else {
        driftDataFrame[elementIndex] <- rep(1, parameters$sampleNumber)
      }
    }
    
    return(list(signal = driftDataFrame, RSD = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber)))
  })
  
  process$driftAndBlankCorrectedMatrix <- reactive({
    
    signalDriftCorrectedRatio <- process$ratio_cor_b()[[1]] / process$driftFactorMatrix()[["signal"]]
    signalDriftCorrectedRSD <- process$ratio_cor_b()[[2]]
    return(list(signal = signalDriftCorrectedRatio, RSD = signalDriftCorrectedRSD))
  })
  
  process$calibrationCoefficientMatrix <- reactive({
    
    calibrationLinearRegressionSlope <- c()
    elementsWithCalibrationIssues <- c()
    
    for (elementIndex in 1:parameters$analyteNumber){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      calibrationSignalUncertaintyConcentration <- getCalibrationData(elementFullName = elementFullName, signal=process$ratio_cor_b(),
                                                                      stdIdentificationColumn=parameters[["categoricalDataAndTime"]][ , "Level"], stdDataFrame = extracted$standard)
      
      calibrationModel <- lm(Concentration ~ 0+Signal, data=as.data.frame(calibrationSignalUncertaintyConcentration))
      
      if (dim(summary(calibrationModel)$coefficients)[1] == 0) {
        calibrationLinearRegressionSlope <- c(calibrationLinearRegressionSlope, NA)
        elementsWithCalibrationIssues <- c(elementsWithCalibrationIssues, elementFullName)
      }
      else
      {
        calibrationLinearRegressionSlope <- c(calibrationLinearRegressionSlope, summary(calibrationModel)$coefficients[1,1])
      }
      
    }
    
    if (!is.null(elementsWithCalibrationIssues)) {
      shinyalert("Calibration alert", paste("Impossible to compute calibration for", paste(elementsWithCalibrationIssues, collapse = ' '), sep = " "), type = "error")
    }
    
    return(calibrationLinearRegressionSlope)
  })
  
  process$concentration <- reactive({
    
    concentration <- t(t(process$driftAndBlankCorrectedMatrix()[["signal"]])*process$calibrationCoefficientMatrix())
    concentration[concentration < 0] <- "<blk"
    concentrationRSD <- process$driftAndBlankCorrectedMatrix()[["RSD"]]
    concentrationRSD[concentration < 0] <- "N/A"
    
    return(list(signal = concentration, RSD = concentrationRSD))
  })
  
  downloadTabData <- reactive({
    
    selectedIndex <- rowIndexInMain$custom[[input$viewConcentrationIndex]]
    
    downloadTabDataSignal <- switch(input$dataReductionStateSelection,
                                     "Counts per second" = process$analyteCountsPerSecond(),
                                     "Internal Standard matrix (adapted for analytes)" = internalStandardMatrixAdaptedToAnalytes()[["signal"]],
                                     "Internal Standard ratio" = process$ratio()[["signal"]],
                                     "Blank matrix" = process$blk_ratio[["signal"]],
                                     "Blank corrected signal (CPS or ratio)" = process$ratio_cor_b()[["signal"]],
                                     "Drift matrix" = process$driftFactorMatrix()[["signal"]],
                                     "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedMatrix()[["signal"]],
                                     "Concentration" = process$concentration()[["signal"]]
    )
    
    downloadTabDataRsd <- switch(input$dataReductionStateSelection,
                                  "Counts per second" = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                  "Internal Standard matrix (adapted for analytes)" = internalStandardMatrixAdaptedToAnalytes()[["RSD"]],
                                  "Internal Standard ratio" = process$ratio()[["RSD"]],
                                  "Blank matrix" = process$blk_ratio[["RSD"]],
                                  "Blank corrected signal (CPS or ratio)" = process$ratio_cor_b()[["RSD"]],
                                  "Drift matrix" = process$driftFactorMatrix()[["RSD"]],
                                  "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedMatrix()[["RSD"]],
                                  "Concentration" = process$concentration()[["RSD"]]
    )
    
    downloadTabDataSignalFormated <- setMatrixFormat(downloadTabDataSignal, columnNamesToAdd = "Sample Name", parameters = parameters)
    downloadTabDataRsdFormated <- setMatrixFormat(downloadTabDataRsd, columnNamesToAdd = "Sample Name", parameters = parameters)
    mergedSignalAndRsdFormated <- setMatrixFormat(list(signal = downloadTabDataSignal, RSD = downloadTabDataRsd),
                                          columnNamesToAdd = "Sample Name", parameters = parameters)
    
    
    customData <- switch(input$viewCustomDataSwitch,
                         "1" = downloadTabDataSignalFormated,
                         "2" = downloadTabDataRsdFormated,
                         "3" = mergedSignalAndRsdFormated
    )
    
    return(customData[selectedIndex, ])
    
  })
  
# File import and viewing -------------------------------------------------
  
  #Text display of imported file
  output$mainFileAssignmentText <- renderText({
    mainFile = uploadedFile$main
    createColoredTextForBooleanViewing(!is.null(mainFile), stateText = "Raw file name: ", invalidStateText = "Unassigned", validStateText = mainFile$name)
  })
  output$standardFileAssignmentText <- renderText({
    standardFile = uploadedFile$standard
    createColoredTextForBooleanViewing(!is.null(standardFile), stateText = "Standard file name: ", invalidStateText = "Unassigned", validStateText = standardFile$name)
  })
  output$internalStandardFileAssignmentText <- renderText({
    internalStandardFile = uploadedFile$internalStandard
    createColoredTextForBooleanViewing(!is.null(internalStandardFile), stateText = "ISTD file name: ", invalidStateText = "Unassigned", validStateText = internalStandardFile$name)
  })
  
  #Table rendering of current uploaded file
  output$uploadedFilePreviewTable <- renderTable({
    req(input$uploadedFile)
    tablePreview = read.table(input$uploadedFile$datapath, sep =";", header = FALSE)
    rowNumber <- min(nrow(tablePreview), input$fileUpload_nrow) 
    columnNumber <- min(ncol(tablePreview), input$fileUpload_ncolumn) 
    tablePreview[1:rowNumber, 1:columnNumber]
  })
  
  #Buttons to set raw, std and ISTD files
  observeEvent(input$setAsMainFile, {
    req(input$uploadedFile)
    uploadedFile$main <- input$uploadedFile
  })
  observeEvent(input$setAsStandardFile, {
    req(input$uploadedFile)
    uploadedFile$standard <- input$uploadedFile
  })
  observeEvent(input$setAsInternalStandardFile, {
    req(input$uploadedFile)
    uploadedFile$internalStandard <- input$uploadedFile
  })
  
  #Button to download the list of analyte and ISTD
  output$downloadISTDTemplate <- downloadHandler(
    filename = "ISTD_Template.csv",
    content = function(file) {
      
      req(isExtractionReady() & !is.null(!is.null(extracted$internalStandard)))
      
      write.csv(createISTDtemplate((uploadedFile$main)$datapath),
                file, sep=";", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
      })

  
  isExtractionReady <- reactive({
    !is.null(uploadedFile$main) & !is.null(uploadedFile$standard)
    })
  
  #Button to extract important information and signal of the dataframe
  observeEvent(input$extract, {
    
    req(isExtractionReady())
    
    mainFileDatapath = uploadedFile$main$datapath
    standardFileDatapath = uploadedFile$standard$datapath
    internalStandardFileDatapath = uploadedFile$internalStandard$datapath
    
    extracted$firstRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep=';')
    extracted$firstRowOfMain <- fillEmptyStrings(extracted$firstRowOfMain)
    
    extracted$secondRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep=';', skip=1)
    
    extracted$main <- read.table(mainFileDatapath, skip = 2, header = FALSE, sep=';', stringsAsFactors=FALSE)
    colnames(extracted$main) <- extracted$firstRowOfMain
    
    extracted$standard <- read.table(standardFileDatapath, header = TRUE, sep=';', stringsAsFactors=FALSE)
    extracted$standard <- removeDuplicateLines(extracted$standard)
    row.names(extracted$standard) <- extracted$standard[,1]
    extracted$standard <- extracted$standard[ , -1]
    
    # columnIndexInMain$analytesCountsPerSecond <- extracted$secondRowOfMain == "CPS" & !grepl("ISTD", extracted$firstRowOfMain)
    # columnIndexInMain$analytesCountsPerSecondRelativeStandardDeviation <- extracted$secondRowOfMain == "CPS RSD" & !grepl("ISTD", extracted$firstRowOfMain)
    # columnIndexInMain$internalStandardCountsPerSecond <- extracted$secondRowOfMain == "CPS" & grepl("ISTD", extracted$firstRowOfMain)
    # columnIndexInMain$internalStandardCountsPerSecondRelativeStandardDeviation <- extracted$secondRowOfMain == "CPS RSD" & grepl("ISTD", extracted$firstRowOfMain)
    # columnIndexInMain$numericalColumns <- min(which(columnIndexInMain$analytesCountsPerSecond)):length(extracted$main)
    
    parameters$sampleNumber <- nrow(extracted$main)
    parameters$analyteNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & !grepl("ISTD", extracted$firstRowOfMain)]
    parameters$analyteNumber <- length(parameters$analyteNames)
    parameters$internalStandardNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & grepl("ISTD", extracted$firstRowOfMain)]
    parameters$internalStandardNumber <- length(parameters$internalStandardNames)

    sampleTime <- as.POSIXct(extracted$main[ , which(extracted$secondRowOfMain == "Acq. Date-Time")], format="%d/%m/%Y %H:%M")
    parameters[["categoricalDataAndTime"]] <- data.frame(sampleTime,
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Sample Name")],
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Type")],
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Level")],
                                                    sampleTime - sampleTime[1],
                                                    stringsAsFactors = FALSE
    )
    colnames(parameters[["categoricalDataAndTime"]]) <- c("Time", "Sample Name", "Type", "Level", "Delta Time")
    
    numericalColumns <- seq(from = max(which(extracted$firstRowOfMain == "Sample")) + 1, to = ncol(extracted$main))
    extracted$main[ , numericalColumns] <- sapply(extracted$main[ , numericalColumns], as.numeric)
    
    if (!is.null(internalStandardFileDatapath)) 
    {
      extracted$internalStandard  <- read.table(internalStandardFileDatapath, header = TRUE, sep= ';', stringsAsFactors = FALSE)
      process$blk_ratio <- propagateUncertainty(a = list(signal = process$analyteCountsPerSecond(), RSD = process$analyteCountsPerSecondRelativeStandardDeviation()), b = internalStandardMatrixAdaptedToAnalytes(), operation="division")
    }
    else 
    {
      process$blk_ratio <- list(signal=process$analyteCountsPerSecond(), RSD = process$analyteCountsPerSecondRelativeStandardDeviation())
    }
    
    rowIndexInMain$custom[["All"]] <- which(rep(x = TRUE, parameters$sampleNumber))
    parameters$driftCorrectedElements <- rep(FALSE, parameters$analyteNumber)
    
    applicationState$isExtractionSuccessful <- TRUE
  })
  
  observeEvent(applicationState$isExtractionSuccessful, {
    if (applicationState$isExtractionSuccessful == TRUE) {
      shinyjs::enable(selector = '.navbar-nav a[data-value="Index creation"')
      shinyjs::enable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
      shinyjs::enable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
      shinyjs::enable(selector = '.navbar-nav a[data-value="Process"')
    }
    else if (applicationState$isExtractionSuccessful == FALSE) {
      shinyjs::disable(selector = '.navbar-nav a[data-value="Index creation"')
      shinyjs::disable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
      shinyjs::disable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
      shinyjs::disable(selector = '.navbar-nav a[data-value="Process"')
    }
  })
  
  #Here we render warnings texts to help the user
  output$extractionReadyText <- renderText({
    createColoredTextForBooleanViewing(isExtractionReady(), stateText = "Data extraction ", invalidStateText = "impossible", validStateText = "ready")
  })
  
# Index creation ----------------------------------------------------------

  output$indexTable <- DT::renderDT(datatable({

    if (is.null(process$analyteCountsPerSecond())){return()}
    
    searchWhere = input$searchIndexwhere
    firstColumnName = ""
    
    if (searchWhere == 'smp'){
      headerRows = parameters[["categoricalDataAndTime"]][ , "Sample Name"]
      firstColumnName = "Sample Names"}
    else if (searchWhere == 'lvl'){
      headerRows = parameters[["categoricalDataAndTime"]][ , "Level"]
      firstColumnName = "Levels"}
    else if (searchWhere == 'type'){
      headerRows = typeColumn()
      firstColumnName = "Type"}
    
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchType == 'ematch'){searchWhat = paste("^", searchWhat, "$", sep="")}
    
    if (searchWhat == ""){
      rowIndexInMain$temp <- c(1:parameters$sampleNumber)
    }
    else{
      rowIndexInMain$temp <- grep(searchWhat, headerRows)
    }
    
    if (displayWhat == "ISTD" & !is.null(process$internalStandardCountsPerSecond())){
      indexTable = cbind(headerRows[rowIndexInMain$temp], process$internalStandardCountsPerSecond()[rowIndexInMain$temp,])
      colnames(indexTable) <- c(firstColumnName, parameters$internalStandardNames)
    }
    else if (displayWhat == "analytes" & !is.null(process$analyteCountsPerSecond())){
      indexTable = cbind(headerRows[rowIndexInMain$temp], process$analyteCountsPerSecond()[rowIndexInMain$temp,])
      colnames(indexTable) <- c(firstColumnName, parameters$analyteNames)
    }
    else {return()}
    
    indexTable
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  observeEvent(input$searchIndexwhat, {
    updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
  })
  
  observeEvent(input$searchIndexCreate, {
    indexName = input$searchIndexName
    customNumIndex = c(1:parameters$sampleNumber)[rowIndexInMain$temp][sort(input$indexTable_rows_selected)]
    if (is.null(rowIndexInMain$custom)) {
      rowIndexInMain$custom <- list()
    }
    rowIndexInMain$custom[[input$searchIndexName]] <- customNumIndex
  })
  
  indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
  
  observeEvent(input$indexSelectAll, {
    if (is.null(sort(input$indexTable_rows_selected))){
      DT::selectRows(indexTableProxy, c(1:length(rowIndexInMain$temp)))
    }
    else if (sort(input$indexTable_rows_selected) == c(1:length(rowIndexInMain$temp))){
      DT::selectRows(indexTableProxy, NULL)
    }
    else{
      DT::selectRows(indexTableProxy, c(1:length(rowIndexInMain$temp)))
    }
  })
  
  observeEvent(rowIndexInMain$custom, {
    if (is.null(rowIndexInMain$custom)){return()}
    updateSelectInput(session,"indexBlkchoiceWhat", label  = "Replace what:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"indexBlkchoiceIn", label  = "Replace in:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"selectDriftIndex", label  = "Define drift index:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(rowIndexInMain$custom)),"All")
  })

# Blank verif/process -----------------------------------------------------

  #Render ISTD table if all conditions are met
  output$blkTable <- DT::renderDT(datatable({
    if (is.null(process$analyteCountsPerSecond())){return()}
    blkMode = input$blkInteractionMode
    
    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]
    
    df_view <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], process$blk_ratio[[1]])
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    
    df_process <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], liveReplaceBlkTable()[[1]][, lc:uc, drop = FALSE])
    names(df_process) <- c("Sample Name", names(df_process)[2:length(df_process)])
    df_process <- format(df_process, digits = 3, scientific=T)

    if (blkMode == "view") {df_view}
    else if (blkMode == "process") {df_process}
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  #Defines a proxy for changing selections in the table
  blkTableProxy <- DT::dataTableProxy("blkTable", session = session)
  
  
  observeEvent(input$blkColSlider, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, rowIndexInMain$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceWhat, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, rowIndexInMain$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceIn, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, rowIndexInMain$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInteractionMode, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, rowIndexInMain$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInterpolationMethod, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, rowIndexInMain$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setBlkInterpolationMethod, {
    repIndex = rowIndexInMain$custom[[input$indexBlkchoiceWhat]]
    if (is.null(liveReplaceBlkTable())| is.null(repIndex)) {return()}

    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]

    process$blk_ratio[[1]][repIndex,lc:uc] <- liveReplaceBlkTable()[[1]][repIndex, lc:uc, drop = FALSE]
    process$blk_ratio[[2]][repIndex,lc:uc] <- liveReplaceBlkTable()[[2]][repIndex, lc:uc, drop = FALSE]
  })
  
  #Updates the slider for row and col selections when modifications are made in process$internalStandardCountsPerSecond(), i.e. when the extract button is hit
  observeEvent(process$analyteCountsPerSecond(), {
    if (is.null(process$analyteCountsPerSecond())){return()}
    updateSliderInput(session,"blkRowSlider", max=parameters$sampleNumber, value = c(1,parameters$sampleNumber))
    updateSliderInput(session,"blkColSlider", max=parameters$analyteNumber, value = c(1,parameters$analyteNumber))
  })
  
  observe({
    if(input$blkInteractionMode == "process") {
      shinyjs::enable("setBlkInterpolationMethod")
      shinyjs::enable("blkColSlider")
      shinyjs::enable("blkInterpolationMethod")
      shinyjs::enable("indexBlkchoiceIn")
    }
    else if(input$blkInteractionMode == "view") {
      shinyjs::disable("setBlkInterpolationMethod")
      shinyjs::disable("blkColSlider")
      shinyjs::disable("blkInterpolationMethod")
      shinyjs::disable("indexBlkchoiceIn")
    }
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$searchBlkSelect, {
    DT::selectRows(liveReplaceBlkTableProxy, NULL)
  })

# Drift verif/process -----------------------------------------------------

  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      rowIndexInMain$drift <- rowIndexInMain$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- DT::renderDT(datatable({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(rowIndexInMain$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- process$ratio_cor_b()[[1]]
    if (input$e_ind_drift < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$e_ind_drift)%%5 == 0))
    }
    uc <- min(lc + 4, parameters$analyteNumber)
    driftTable <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"][rowIndexInMain$drift],driftSignal[rowIndexInMain$drift,lc:uc])
    names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
    driftTable
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  
  observeEvent(parameters$analyteNames, {
    if (is.null(parameters$analyteNames)){return()}
    updateSelectInput(session,"e_drift",choices=parameters$analyteNames,selected=parameters$analyteNames[1])
    updateNumericInput(session,"e_ind_drift", "Element number", 1, min = 1, max = length(parameters$analyteNames))
  })
  
  observeEvent(input$changeElementDisplayedForDrift, {
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,parameters$analyteNames,fixed=TRUE))
  })
  
  observeEvent(input$e_ind_drift, {
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=parameters$analyteNames[input$e_ind_drift])
  })
  
  output$driftPlot <- renderPlot({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(rowIndexInMain$drift)){return()}
    elementName = parameters$analyteNames[input$e_ind_drift]
    driftTime = parameters[["categoricalDataAndTime"]][ , "Delta Time"][rowIndexInMain$drift]
    driftValue = process$ratio_cor_b()[[1]][rowIndexInMain$drift,input$e_ind_drift]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = elementName, pch=21, bg="lightblue")

    if (parameters$driftCorrectedElements[input$e_ind_drift] == TRUE){
      time_0 = driftTime[1]
      time_f = driftTime[length(driftTime)]
      timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
      
      dt = as.numeric(driftTime)
      driftModel <- lm(driftValue ~ poly(dt, degree=2, raw=TRUE))
      
      driftPredict=predict(driftModel, newdata = data.frame(dt = timeInterval))
      
      lines(x=timeInterval, y=driftPredict, col="red")
    }
    arrows(parameters[["categoricalDataAndTime"]][ , "Delta Time"][rowIndexInMain$drift], yminus(), parameters[["categoricalDataAndTime"]][ , "Delta Time"][rowIndexInMain$drift], yplus(), length=0.05, angle=90, code=3)
  })
  
  observeEvent(input$setDriftCorrection, {
    if (is.null(parameters$driftCorrectedElements)){return()}
    parameters$driftCorrectedElements[input$e_ind_drift] <- !parameters$driftCorrectedElements[input$e_ind_drift]
  })
  
  observeEvent(input$pressedKey, {
    if (is.null(parameters$driftCorrectedElements) || input$pressedKeyId != C_LETTER_KEYCODE || input$tagId != "e_ind_drift"){return()}
    parameters$driftCorrectedElements[input$e_ind_drift] <- !parameters$driftCorrectedElements[input$e_ind_drift]
  })

# Process -----------------------------------------------------------------
  
  output$customDataTable <- DT::renderDT(datatable({
    
    downloadTabData()
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  

# Downloadable data -------------------------------------------------------
 
  output$downloadCustomTable <- downloadHandler(
    filename = paste0("data_", input$viewConcentrationIndex, ".csv"),
    content = function(file) {
      
      write.csv(downloadTabData(), file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
  )
}
