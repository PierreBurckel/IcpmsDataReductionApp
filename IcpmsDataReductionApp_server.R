
C_LETTER_KEYCODE <- 67

ICPMS_server <- function(input, output, session) {
  
  printLogJs <- function(x, ...) {
    
    logjs(x)
    
    T
  }
  
  addHandler(printLogJs)
  
# Code to be executed at start --------------------------------------------
  
  shinyjs::disable(selector = '.navbar-nav a[data-value="Index creation"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Process"')

# Declaration of reactive values ------------------------------------------
  
  uploadedFile <- reactiveValues()
  extracted <- reactiveValues()
  rowIndexInMain <- reactiveValues()
  process <- reactiveValues()
  modifiers <- reactiveValues(blank = list())
  parameters <- reactiveValues()
  applicationState <- reactiveValues(isExtractionSuccessful = FALSE)

# Reactive expressions for data reduction ---------------------------------
  
  process$countsPerSecond <- reactive({
    countsPerSecond <- extracted$main[ , extracted$secondRowOfMain == "CPS", drop=FALSE]
    numericalCountsPerSecond <- apply(countsPerSecond, c(1,2), is.numeric)
    countsPerSecond[!numericalCountsPerSecond] <- NA
    as.matrix(countsPerSecond)
  })
  
  process$analyteCountsPerSecond <- reactive({
    process$countsPerSecond()[ , parameters$analyteNames, drop=FALSE]
  })
  
  process$internalStandardCountsPerSecond <- reactive({
    process$countsPerSecond()[ , parameters$internalStandardNames, drop=FALSE]
  })
  
  process$relativeStandardDeviation <- reactive({
    relativeStandardDeviation <- extracted$main[ , extracted$secondRowOfMain == "CPS RSD", drop=FALSE]
    numericalRelativeStandardDeviation <- apply(relativeStandardDeviation, c(1,2), is.numeric)
    relativeStandardDeviation[!numericalRelativeStandardDeviation] <- NA
    as.matrix(relativeStandardDeviation)
  })
  
  process$analyteCountsPerSecondRelativeStandardDeviation <- reactive({
    process$relativeStandardDeviation()[ , parameters$analyteNames, drop=FALSE]
  })
  
  process$internalStandardCountsPerSecondRelativeStandardDeviation <- reactive({
    process$relativeStandardDeviation()[ , parameters$internalStandardNames, drop=FALSE]
  })
  
  process$analyteCountsPerSecondEudc <- reactive({
    EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                        estimatedData = process$analyteCountsPerSecond(),
                                        uncertaintyData = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                        uncertaintyType = "rsd")
  })
  
  process$internalStandardCountsPerSecondEudc <- reactive({
    EstimationUncertaintyDataCouple$new(elementFullNames = parameters$internalStandardNames,
                                        estimatedData = process$internalStandardCountsPerSecond(),
                                        uncertaintyData = process$internalStandardCountsPerSecondRelativeStandardDeviation(),
                                        uncertaintyType = "rsd")
  })
  
  process$internalStandardEudcAdaptedToAnalytes <- reactive({
    if (!is.null(extracted$internalStandardToAnalyteAssignment)) 
    {
      return(createinternalStandardEudcAdaptedToAnalytes(internalStandardToAnalyteAssignmentDataframe = extracted$internalStandardToAnalyteAssignment,
                                                         internalStandardCountsPerSecondEudc = process$internalStandardCountsPerSecondEudc(),
                                                         parameters = parameters))
    }
    else
    {
      return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                 estimatedData = matrix(1, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                                 uncertaintyData = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                                 uncertaintyType = "sd"))
    }
  })
  
  # process$interferenceCorrectedCountsPerSecondEudc <- reactive({
  #   for (interference in parameters$interferenceParametersList) {
  #     interferenceCorrectedCountsPerSecondEudc <- correctInterferences(eudcToCorrect = process$analyteCountsPerSecondEudc(), interferenceParameters = interference)
  #   }
  #   return(interferenceCorrectedCountsPerSecondEudc)
  # })
  
  process$interferenceCorrectedCountsPerSecond <- reactive({
    process$analyteCountsPerSecondEudc()$subtract(process$interferenceCountsPerSeconds)
  })
  
  process$analyteToIstdRatio <- reactive({
    process$interferenceCorrectedCountsPerSecond()$divideBy(process$internalStandardEudcAdaptedToAnalytes())
    })
  
  process$analyteToIstdBlankRatio <- reactive({
    if (input$enableBlankCorrection == TRUE)
    {
      process$analyteToIstdRatio() %>% applyModifierToEudc(modifiers$blank)
    }
    else
    {
      EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                          estimatedData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                          uncertaintyData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                          uncertaintyType = "sd")
    }
  })
  
  process$analyteToIstdRatioBlankCorrected <- reactive({
    process$analyteToIstdRatio()$subtract(process$analyteToIstdBlankRatio())
  })
  
  parameters$deltaTime <- reactive({
    deltaTime <- as.numeric(parameters[["categoricalDataAndTime"]][ , "Time"] - 
                            parameters[["categoricalDataAndTime"]][1, "Time"])
    return(deltaTime)
  })
  
  parameters$listOfElementSpecificDriftIndex <- reactive({
    
    req(parameters$analyteNames)
    
    listOfElementSpecificDriftIndex <- list()
    
    for (elementFullName in parameters$analyteNames) {
      listOfElementSpecificDriftIndex[[elementFullName]] <- getElementSpecificDriftIndex(elementFullName = elementFullName,
                                                                                         standardDataFrame = extracted$standard, 
                                                                                         standardIdentificationColumn = parameters[["categoricalDataAndTime"]][ , "Level"],
                                                                                         driftIndex = rowIndexInMain$drift)
    }
    
    return(listOfElementSpecificDriftIndex)
  })
  
  process$elementSpecificDriftModels <- reactive({
    
    listOfDriftModels <- list()
  
    for (elementIndex in seq(parameters$analyteNumber)){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      elementSpecificDriftIndex <- parameters$listOfElementSpecificDriftIndex()[[elementFullName]]
      
      driftSignalAndDeltaTime <- cbind(Signal=process$analyteToIstdRatioBlankCorrected()$getEstimation()[elementSpecificDriftIndex,elementIndex],
                                       dt=parameters$deltaTime()[elementSpecificDriftIndex])  
      
      if (length(elementSpecificDriftIndex) < 3){
        
        driftModel = NA
      }
      else{
        driftModel <- lm(Signal ~ poly(dt, degree=2, raw=TRUE), data = as.data.frame(driftSignalAndDeltaTime))
      }
      
      listOfDriftModels[[elementFullName]] <- driftModel
    }
    
    return(listOfDriftModels)
  })
    
  process$driftFactorEudc <- reactive({
    for (elementIndex in seq(parameters$analyteNumber)){
        
      elementFullName <- parameters$analyteNames[elementIndex]
      elementSpecificDriftIndex <- parameters$listOfElementSpecificDriftIndex()[[elementFullName]]
      
      if (all(is.na(process$elementSpecificDriftModels()[[elementFullName]]))){
        driftPredict = rep(NA, parameters$sampleNumber)
      }
      else{
        driftModel <- process$elementSpecificDriftModels()[[elementFullName]]
        driftPredict=predict(driftModel, newdata = data.frame(dt=parameters$deltaTime()))
      }
      
      if (elementIndex == 1){
        driftMatrix <- as.matrix(driftPredict)
      }
      else{
        driftMatrix <- cbind(driftMatrix, driftPredict)
      }
  
      if (parameters$driftCorrectedElements[elementIndex] == TRUE){
        driftMatrix[ , elementIndex] <- driftMatrix[ , elementIndex] / driftMatrix[min(elementSpecificDriftIndex), elementIndex]
        driftMatrix[is.na(driftMatrix[ , elementIndex]), elementIndex] <- 1
        driftMatrix[seq(min(elementSpecificDriftIndex)), elementIndex] <- 1
      } 
      else {
        driftMatrix[ , elementIndex] <- rep(1, parameters$sampleNumber)
      }
    }
    return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                               estimatedData = driftMatrix,
                                               uncertaintyData = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                               uncertaintyType = "sd"))
  })
  
  process$driftAndBlankCorrectedEudc <- reactive({
    driftAndBlankCorrectedEudc <- process$analyteToIstdRatioBlankCorrected()$divideBy(process$driftFactorEudc())
    return(driftAndBlankCorrectedEudc)
  })
  
  process$calibrationCoefficientMatrix <- reactive({
    
    calibrationLinearRegressionSlope <- c()
    elementsWithCalibrationIssues <- c()
    
    for (elementIndex in 1:parameters$analyteNumber){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      calibrationSignalUncertaintyConcentration <- getCalibrationData(elementFullName = elementFullName, signalMatrix = process$analyteToIstdRatioBlankCorrected(),
                                                                      standardIdentificationColumn = parameters[["categoricalDataAndTime"]][ , "Level"], standardDataFrame = extracted$standard)
      
      if (!is.null(calibrationSignalUncertaintyConcentration)) {
        calibrationModel <- lm(Concentration ~ 0+Signal, data=as.data.frame(calibrationSignalUncertaintyConcentration))
      }
      
      if (is.null(calibrationSignalUncertaintyConcentration) || dim(summary(calibrationModel)$coefficients)[1] == 0) {
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
    
    concentration <- t(t(process$driftAndBlankCorrectedEudc()$getEstimation())*process$calibrationCoefficientMatrix())
    concentrationRSD <- process$driftAndBlankCorrectedEudc()$getRsd()
    
    return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                               estimatedData = concentration,
                                               uncertaintyData = concentrationRSD,
                                               uncertaintyType = "rsd"))
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
    tablePreview = read.table(input$uploadedFile$datapath, sep =input$csvDelimitation, header = FALSE)
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
      
      req(isExtractionReady() & !is.null(!is.null(extracted$internalStandardToAnalyteAssignment)))
      
      write.table(createISTDtemplate(dataFileName = (uploadedFile$main)$datapath,
                                     sep = input$csvDelimitation),
                file, sep = input$csvDelimitation, quote = FALSE,
                row.names = FALSE, col.names = TRUE)
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
    
    extracted$firstRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation)
    extracted$firstRowOfMain <- fillEmptyStrings(extracted$firstRowOfMain)
    
    extracted$secondRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation, skip=1)
    
    requiredStringsInFirstRowOfMain <- "Sample"
    requiredStringsInSecondRowOfMain <- c("Acq. Date-Time", "Sample Name", "Type", "Level", "CPS", "CPS RSD")
    
    if (!all(requiredStringsInFirstRowOfMain %in% extracted$firstRowOfMain) || !all(requiredStringsInSecondRowOfMain %in% extracted$secondRowOfMain)) {
      shinyalert("Impossible to extract", paste("Missing either the Sample column in the first row of the main file, or  one of the following columns in the second row: ", paste(requiredStringsInSecondRowOfMain, collapse = ', '), sep = " "), type = "error")
      return(NULL)
    }
    
    extracted$main <- read.table(mainFileDatapath, skip = 2, header = FALSE, sep = input$csvDelimitation, stringsAsFactors=FALSE)
    colnames(extracted$main) <- extracted$firstRowOfMain
    
    extracted$standard <- createStandardDataFrameFromFile(dataPath = standardFileDatapath, sep = input$csvDelimitation)
    if (is.null(extracted$standard)) return(NULL)
    
    parameters$sampleNumber <- nrow(extracted$main)
    parameters$analyteNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & !grepl("ISTD", extracted$firstRowOfMain)]
    parameters$analyteNumber <- length(parameters$analyteNames)
    parameters$internalStandardNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & grepl("ISTD", extracted$firstRowOfMain)]
    parameters$internalStandardNumber <- length(parameters$internalStandardNames)
    
    sampleTime <- as.POSIXct(extracted$main[ , which(extracted$secondRowOfMain == "Acq. Date-Time")], format= input$dateFormat)
    if(any(is.na(sampleTime))) {
      shinyalert("Impossible to extract", paste0("Error, non-valid (NA) dates at sample number(s) ", paste(which(is.na(sampleTime)), collapse = ", "),  ". Check format in csv file and in File upload and parameters tab"), type = "error")
      return(NULL)
    }
    
    parameters[["categoricalDataAndTime"]] <- data.frame(sampleTime,
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Sample Name")],
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Type")],
                                                    extracted$main[ , which(extracted$secondRowOfMain == "Level")],
                                                    stringsAsFactors = FALSE
    )
    colnames(parameters[["categoricalDataAndTime"]]) <- c("Time", "Sample Name", "Type", "Level")
    
    numericalColumns <- seq(from = max(which(extracted$firstRowOfMain == "Sample")) + 1, to = ncol(extracted$main))
    extracted$main[ , numericalColumns] <- sapply(extracted$main[ , numericalColumns], as.numeric)
    
    process$interferenceCountsPerSeconds <- EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                                                    estimatedData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                                                                    uncertaintyData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                                                                    uncertaintyType = "sd")
    
    # if (!is.null(internalStandardFileDatapath)) 
    # {
    #   extracted$internalStandardToAnalyteAssignment  <- read.table(internalStandardFileDatapath, header = TRUE, sep= input$csvDelimitation, stringsAsFactors = FALSE)
    #   process$analyteToIstdBlankRatio <- process$analyteCountsPerSecondEudc()$divideBy(process$internalStandardEudcAdaptedToAnalytes())
    # }
    # else 
    # {
    #   process$analyteToIstdBlankRatio <- process$analyteCountsPerSecondEudc()
    # }
    
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
      headerRows = parameters[["categoricalDataAndTime"]][ , "Type"]
      firstColumnName = "Type"}
    
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchType == 'ematch'){searchWhat = paste("^", searchWhat, "$", sep="")}
    
    if (searchWhat == ""){
      rowIndexInMain$index_rowsMatchingRegularExpression <- seq(parameters$sampleNumber)
    }
    else{
      rowIndexInMain$index_rowsMatchingRegularExpression <- grep(searchWhat, headerRows)
    }
    
    if (displayWhat == "ISTD" & !is.null(process$internalStandardCountsPerSecondEudc()$getEstimation())){
      indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process$internalStandardCountsPerSecondEudc()$getEstimation()[rowIndexInMain$index_rowsMatchingRegularExpression,])
      colnames(indexTable) <- c(firstColumnName, parameters$internalStandardNames)
    }
    else if (displayWhat == "analytes" & !is.null(process$analyteCountsPerSecondEudc()$getEstimation())){
      indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process$analyteCountsPerSecondEudc()$getEstimation()[rowIndexInMain$index_rowsMatchingRegularExpression,])
      colnames(indexTable) <- c(firstColumnName, parameters$analyteNames)
    }
    else {return()}
    
    indexTable
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  indexTab_selectedRows <- reactive({
    sort(input$indexTable_rows_selected)
  })
  
  observeEvent(input$searchIndexwhat, {
    updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
  })
  
  observeEvent(input$searchIndexCreate, {
    indexName = input$searchIndexName
    customNumIndex = c(1:parameters$sampleNumber)[rowIndexInMain$index_rowsMatchingRegularExpression][indexTab_selectedRows()]
    if (is.null(rowIndexInMain$custom)) {
      rowIndexInMain$custom <- list()
    }
    rowIndexInMain$custom[[input$searchIndexName]] <- customNumIndex
  })
  
  indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
  
  observeEvent(input$indexSelectAll, {
    
    if (is.null(indexTab_selectedRows()))
    {
      DT::selectRows(indexTableProxy, seq(rowIndexInMain$index_rowsMatchingRegularExpression))
    }
    else if (identical(indexTab_selectedRows(), seq(rowIndexInMain$index_rowsMatchingRegularExpression)))
    {
      DT::selectRows(indexTableProxy, NULL)
    }
    else
    {
      DT::selectRows(indexTableProxy, seq(rowIndexInMain$index_rowsMatchingRegularExpression))
    }
    
  })
  
  observeEvent(rowIndexInMain$custom, {
    if (is.null(rowIndexInMain$custom)){return()}
    updateSelectInput(session,"sliderInput_InterferenceTab_indexForCorrectionFactorCalculation", label  = "Compute correction factor from:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"sliderInput_BlankTab_rowsToReplace", label  = "Replace what:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"sliderInput_BlankTab_rowsToReplaceFrom", label  = "Replace in:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"selectDriftIndex", label  = "Define drift index:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(rowIndexInMain$custom)),"All")
  })


# Interference correction -------------------------------------------------

  observeEvent(input$actionButton_InterferenceTab_applyCorrection, {
    
    interferenceSignal <- process$interferenceCountsPerSeconds$getEstimation()
    interferenceRsd <- process$interferenceCountsPerSeconds$getRsd()
    analyteSignal <- process$analyteCountsPerSecondEudc()$getEstimation()
    analyteRsd <- process$analyteCountsPerSecondEudc()$getRsd()
    
    indexForCorrectionFactorCalculation <- rowIndexInMain$custom[[input$sliderInput_InterferenceTab_indexForCorrectionFactorCalculation]]
    interferedElement <- input$sliderInput_InterferenceTab_interferedElement
    interferingElement <- input$sliderInput_InterferenceTab_interferingElement
    
    if (length(indexForCorrectionFactorCalculation) > 1) {
      estimatedCorrectionFactors <- analyteSignal[indexForCorrectionFactorCalculation, interferedElement] / analyteSignal[indexForCorrectionFactorCalculation, interferingElement]
      correctionFactorValue <- mean(estimatedCorrectionFactors)
      correctionFactorRsd <- sd(estimatedCorrectionFactors) / correctionFactorValue * 100
    }
    else {
      correctionFactorValue <- analyteSignal[indexForCorrectionFactorCalculation, interferedElement] / analyteSignal[indexForCorrectionFactorCalculation, interferingElement]
      correctionFactorRsd <- sqrt((analyteRsd[indexForCorrectionFactorCalculation, interferedElement]/100)^2 +
                                  (analyteRsd[indexForCorrectionFactorCalculation, interferingElement]/100)^2) * 100
    }
    
    interferenceSignal[ , interferedElement] <- analyteSignal[ , interferingElement] * correctionFactorValue
    interferenceRsd[ , interferedElement] <- sqrt((analyteRsd[ , interferingElement]/100)^2 + (correctionFactorRsd/100)^2) * 100
    
    newInterferenceEudc <- EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                           estimatedData = interferenceSignal,
                                                           uncertaintyData = interferenceRsd,
                                                           uncertaintyType = "rsd")
    
    process$interferenceCountsPerSeconds <- process$interferenceCountsPerSeconds$add(newInterferenceEudc)
  })  
  
# Blank verif/process -----------------------------------------------------

  activeBlankModifier <- reactive({
    list(DataModifier$new(linesToBeReplaced = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]], 
                          columnsToBeReplaced = 1:parameters$analyteNumber,
                          linesUsedForReplacement = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplaceFrom]],
                          columnsUsedForReplacement = 1:parameters$analyteNumber,
                          howToReplace = input$sliderInput_BlankTab_replacementMethod))
  })
  # liveReplaceBlkTable <- reactive({
  #   
  #   process$analyteToIstdRatio()$applyModifications(activeBlankModifier())
  #   
  #   # replaceIndexWhat = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]]
  #   # replaceIndexIn = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplaceFrom]]
  #   # 
  #   # replaceValues(process$analyteToIstdRatio(), input$sliderInput_BlankTab_replacementMethod, replaceIndexWhat, replaceIndexIn, parameters)
  # })
  
  #Render ISTD table if all conditions are met
  output$blankTab_table <- DT::renderDT(datatable({
    
    blank_processOrView = input$sliderInput_BlankTab_processOrView
    
    if (blank_processOrView == "view") 
    {
      modifiedAnalyteToIstdRatio <- process$analyteToIstdBlankRatio()
      blankTab_table <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], modifiedAnalyteToIstdRatio$getEstimation())
      names(blankTab_table) <- c("Sample Name", names(blankTab_table)[2:length(blankTab_table)])
      blankTab_table <- format(blankTab_table, digits = 3, scientific=T)
    }
    if (blank_processOrView == "process") 
    {
      activelyModifiedAnalyteToIstdRatio <- process$analyteToIstdRatio() %>% applyModifierToEudc(activeBlankModifier())
      blankTab_table <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], activelyModifiedAnalyteToIstdRatio$getEstimation())
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
      DT::selectRows(blankTab_tableProxy, rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]])
    }
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$actionButton_BlankTab_replace, {
    blankModifierNumber <- length(modifiers$blank)
    modifiers$blank[blankModifierNumber + 1] <- activeBlankModifier()
    
    # rowReplacementIndex = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]]
    # 
    # req(!is.null(liveReplaceBlkTable()) && !is.null(rowReplacementIndex))
    # 
    # blankRatioSignal <- process$analyteToIstdBlankRatio$getEstimation()
    # blankRatioRsd <- process$analyteToIstdBlankRatio$getRsd()
    # blankRatioSignal[rowReplacementIndex, ] <- liveReplaceBlkTable()$getEstimation()[rowReplacementIndex, , drop = FALSE]
    # blankRatioRsd[rowReplacementIndex, ] <- liveReplaceBlkTable()$getRsd()[rowReplacementIndex, , drop = FALSE]
    # 
    # process$analyteToIstdBlankRatio <- EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
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

# Drift verif/process -----------------------------------------------------

  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      rowIndexInMain$drift <- rowIndexInMain$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- DT::renderDT(datatable({
    if (is.null(process$analyteToIstdRatioBlankCorrected()$getEstimation()) | is.null(rowIndexInMain$drift) | !is.integer(input$driftTab_numericInput_analyteNumber)){return()}
    driftSignal <- process$analyteToIstdRatioBlankCorrected()$getEstimation()
    if (input$driftTab_numericInput_analyteNumber < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$driftTab_numericInput_analyteNumber)%%5 == 0))
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
  
  yplus <- reactive({(process$analyteToIstdRatioBlankCorrected()$getEstimation() + process$analyteToIstdRatioBlankCorrected()$getSd())[rowIndexInMain$drift, input$driftTab_numericInput_analyteNumber]})
  yminus <- reactive({(process$analyteToIstdRatioBlankCorrected()$getEstimation() - process$analyteToIstdRatioBlankCorrected()$getSd())[rowIndexInMain$drift, input$driftTab_numericInput_analyteNumber]})
  
  output$driftPlot <- renderPlot({
    if (is.null(process$analyteToIstdRatioBlankCorrected()$getEstimation()) | is.null(rowIndexInMain$drift)){return()}
    
    elementFullName <- parameters$analyteNames[input$driftTab_numericInput_analyteNumber]
    
    driftTime <- parameters$deltaTime()[rowIndexInMain$drift]
    driftValue <- process$analyteToIstdRatioBlankCorrected()$getEstimation()[rowIndexInMain$drift,input$driftTab_numericInput_analyteNumber]
    driftSd <- process$analyteToIstdRatioBlankCorrected()$getSd()[rowIndexInMain$drift,input$driftTab_numericInput_analyteNumber]
    
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

# Process -----------------------------------------------------------------
  
  downloadTabData <- reactive({
    
    selectedIndex <- rowIndexInMain$custom[[input$viewConcentrationIndex]]
    
    downloadTabDataSignal <- switch(input$dataReductionStateSelection,
                                    "Counts per second" = process$analyteCountsPerSecond(),
                                    "Interference counts per second" = process$interferenceCountsPerSeconds$getEstimation(),
                                    "Interference corrected counts per second" = process$interferenceCorrectedCountsPerSecond()$getEstimation(),
                                    "Internal Standard matrix (adapted for analytes)" = process$internalStandardEudcAdaptedToAnalytes()$getEstimation(),
                                    "Internal Standard ratio" = process$analyteToIstdRatio()$getEstimation(),
                                    "Blank matrix" = process$analyteToIstdBlankRatio()$getEstimation(),
                                    "Blank corrected signal (CPS or ratio)" = process$analyteToIstdRatioBlankCorrected()$getEstimation(),
                                    "Drift matrix" = process$driftFactorEudc()$getEstimation(),
                                    "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedEudc()$getEstimation(),
                                    "Concentration" = process$concentration()$getEstimation()
    )
    
    downloadTabDataSignal[downloadTabDataSignal < 0] <- "<blk"
    downloadTabDataSignal[is.na(downloadTabDataSignal)] <- "N/A"
    
    downloadTabDataRsd <- switch(input$dataReductionStateSelection,
                                 "Counts per second" = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                 "Interference counts per second" = process$interferenceCountsPerSeconds$getRsd(),
                                 "Interference corrected counts per second" = process$interferenceCorrectedCountsPerSecond()$getRsd(),
                                 "Internal Standard matrix (adapted for analytes)" = process$internalStandardEudcAdaptedToAnalytes()$getRsd(),
                                 "Internal Standard ratio" = process$analyteToIstdRatio()$getRsd(),
                                 "Blank matrix" = process$analyteToIstdBlankRatio()$getRsd(),
                                 "Blank corrected signal (CPS or ratio)" = process$analyteToIstdRatioBlankCorrected()$getRsd(),
                                 "Drift matrix" = process$driftFactorEudc()$getRsd(),
                                 "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedEudc()$getRsd(),
                                 "Concentration" = process$concentration()$getRsd()
    )
    
    downloadTabDataRsd[downloadTabDataSignal < 0] <- "N/A"
    downloadTabDataRsd[is.na(downloadTabDataRsd)] <- "N/A"
    
    colnames(downloadTabDataSignal) <- parameters[["analyteNames"]]
    colnames(downloadTabDataRsd) <- parameters[["analyteNames"]]
    
    mergedSignalAndRsd <- mergeMatrixes(downloadTabDataSignal, downloadTabDataRsd)
    mergedSignalAndRsdColumnNames <- mergeMatrixes(matrix(parameters[["analyteNames"]], nrow = 1, ncol = parameters[["analyteNumber"]]),
                                                           matrix(paste0(parameters[["analyteNames"]], " RSD (%)"), nrow = 1, ncol = parameters[["analyteNumber"]]))
    colnames(mergedSignalAndRsd) <- mergedSignalAndRsdColumnNames
    
    
    customData <- switch(input$viewCustomDataSwitch,
                         "1" = downloadTabDataSignal,
                         "2" = downloadTabDataRsd,
                         "3" = mergedSignalAndRsd
    )
    
    customDataWithHeader <- cbind(parameters[["categoricalDataAndTime"]][, "Sample Name", drop = FALSE], customData)
    
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
      write.table(downloadTabData(), file, sep=input$csvDelimitation, quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  )
}
