
C_LETTER_KEYCODE <- 67

ICPMS_server <- function(input, output, session) {

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
  parameters <- reactiveValues()
  applicationState <- reactiveValues(isExtractionSuccessful = FALSE)

# Reactive expressions for data reduction ---------------------------------
  
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
  
  process$internalStandardMatrixAdaptedToAnalytes <- reactive({
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
  
  process$ratio <- reactive({
    propagateUncertainty(a = list(signal=process$analyteCountsPerSecond(), RSD=process$analyteCountsPerSecondRelativeStandardDeviation()),
                         b = process$internalStandardMatrixAdaptedToAnalytes(),
                         operation="division")
    })
  
  process$ratio_cor_b <- reactive({
      propagateUncertainty(a = process$ratio(),
                           b = process$blk_ratio,
                           operation="substraction")
  })
  
  parameters$deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- reactive({
    
    listOfdeltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- list()
    
    for (elementIndex in 1:parameters$analyteNumber){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      driftIndexAfterFirstStandard <- getElementSpecificDriftIndex(elementFullName = elementFullName, stdDataFrame = extracted$standard, 
                                                                   stdIdentificationColumn=parameters[["categoricalDataAndTime"]][ , "Level"], driftIndex = rowIndexInMain$drift)
      
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- as.numeric(parameters[["categoricalDataAndTime"]][ , "Time"] - parameters[["categoricalDataAndTime"]][ , "Time"][driftIndexAfterFirstStandard][1])
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[deltaTimeWithFirstDriftAfterFirstStandardAsOrigin < 0] <- 0
      
      listOfdeltaTimeWithFirstDriftAfterFirstStandardAsOrigin[[elementFullName]] <- deltaTimeWithFirstDriftAfterFirstStandardAsOrigin
    }
    
    return(listOfdeltaTimeWithFirstDriftAfterFirstStandardAsOrigin)
  })
  
  process$elementSpecificDriftModels <- reactive({
    
    listOfDriftModels <- list()
  
    for (elementIndex in 1:parameters$analyteNumber){
      
      elementFullName <- parameters$analyteNames[elementIndex]
      
      driftIndexAfterFirstStandard <- getElementSpecificDriftIndex(elementFullName = elementFullName, stdDataFrame = extracted$standard, 
                                                           stdIdentificationColumn=parameters[["categoricalDataAndTime"]][ , "Level"], driftIndex = rowIndexInMain$drift)
      
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- parameters$deltaTimeWithFirstDriftAfterFirstStandardAsOrigin()[[elementFullName]]
      
      driftSignalAndDeltaTime <- cbind(Signal=process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex],
                                       dt=deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[driftIndexAfterFirstStandard])  
      
      if (all(is.na(process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex]))){
        driftModel = NA
      }
      else{
        driftModel <- lm(Signal ~ poly(dt, degree=2, raw=TRUE), data = as.data.frame(driftSignalAndDeltaTime))
      }
      
      listOfDriftModels[[elementFullName]] <- driftModel
    }
    
    return(listOfDriftModels)
  })
    
  process$driftFactorMatrix <- reactive({
      
    for (elementIndex in 1:parameters$analyteNumber){
        
      elementFullName <- parameters$analyteNames[elementIndex]
      
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- parameters$deltaTimeWithFirstDriftAfterFirstStandardAsOrigin()[[elementFullName]]
      
      if (is.na(process$elementSpecificDriftModels()[[elementFullName]])){
        driftPredict = rep(NA, parameters$sampleNumber)
      }
      else{
        driftModel <- process$elementSpecificDriftModels()[[elementFullName]]
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
      
      calibrationSignalUncertaintyConcentration <- getCalibrationData(isotopeName = elementFullName, signalMatrix = process$ratio_cor_b(),
                                                                      standardIdentificationColumn = parameters[["categoricalDataAndTime"]][ , "Level"], standardDataMatrix = extracted$standard)
      
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
    
    concentration <- t(t(process$driftAndBlankCorrectedMatrix()[["signal"]])*process$calibrationCoefficientMatrix())
    concentration[concentration < 0] <- "<blk"
    concentrationRSD <- process$driftAndBlankCorrectedMatrix()[["RSD"]]
    concentrationRSD[concentration < 0] <- "N/A"
    
    return(list(signal = concentration, RSD = concentrationRSD))
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
    
    requiredStringsInFirstRowOfMain <- "Sample"
    requiredStringsInSecondRowOfMain <- c("Acq. Date-Time", "Sample Name", "Type", "Level", "CPS", "CPS RSD")
    
    if (!all(requiredStringsInFirstRowOfMain %in% extracted$firstRowOfMain) || !all(requiredStringsInSecondRowOfMain %in% extracted$secondRowOfMain)) {
      shinyalert("Impossible to extract", paste("Missing either the Sample column in the first row of the main file, or  one of the following columns in the second row: ", paste(requiredStringsInSecondRowOfMain, collapse = ', '), sep = " "), type = "error")
      return(NULL)
    }
    
    extracted$main <- read.table(mainFileDatapath, skip = 2, header = FALSE, sep=';', stringsAsFactors=FALSE)
    colnames(extracted$main) <- extracted$firstRowOfMain
    
    extracted$standard <- read.table(standardFileDatapath, header = FALSE, sep=';', stringsAsFactors=FALSE)
    if (!identical(unique(extracted$standard[1, ]), extracted$standard[1, ])) {
      shinyalert("Impossible to extract", "In the standard file, the standard labels in the first row are not unique", type = "error")
      return(NULL)
    }
    colnames(extracted$standard) <- extracted$standard[1, ]
    extracted$standard <- extracted$standard[-1, ]
    extracted$standard <- removeDuplicateLines(extracted$standard)
    row.names(extracted$standard) <- extracted$standard[ , 1]
    extracted$standard <- extracted$standard[ , -1]
    
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
      process$blk_ratio <- propagateUncertainty(a = list(signal = process$analyteCountsPerSecond(), RSD = process$analyteCountsPerSecondRelativeStandardDeviation()),
                                                b = process$internalStandardMatrixAdaptedToAnalytes(),
                                                operation="division")
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
    
    if (displayWhat == "ISTD" & !is.null(process$internalStandardCountsPerSecond())){
      indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process$internalStandardCountsPerSecond()[rowIndexInMain$index_rowsMatchingRegularExpression,])
      colnames(indexTable) <- c(firstColumnName, parameters$internalStandardNames)
    }
    else if (displayWhat == "analytes" & !is.null(process$analyteCountsPerSecond())){
      indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process$analyteCountsPerSecond()[rowIndexInMain$index_rowsMatchingRegularExpression,])
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
    updateSelectInput(session,"sliderInput_BlankTab_rowsToReplace", label  = "Replace what:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"sliderInput_BlankTab_rowsToReplaceFrom", label  = "Replace in:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"selectDriftIndex", label  = "Define drift index:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(rowIndexInMain$custom)),"All")
  })

# Blank verif/process -----------------------------------------------------

  liveReplaceBlkTable <- reactive({
    
    replaceIndexWhat = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]]
    replaceIndexIn = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplaceFrom]]
    
    replaceValues(process$ratio(), 1:parameters$analyteNumber, input$sliderInput_BlankTab_replacementMethod, replaceIndexWhat, replaceIndexIn)
  })
  
  #Render ISTD table if all conditions are met
  output$blankTab_table <- DT::renderDT(datatable({
    
    blank_processOrView = input$sliderInput_BlankTab_processOrView
    
    if (blank_processOrView == "view") 
    {
      blankTab_table <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], process$blk_ratio[[1]])
      names(blankTab_table) <- c("Sample Name", names(blankTab_table)[2:length(blankTab_table)])
      blankTab_table <- format(blankTab_table, digits = 3, scientific=T)
    }
    if (blank_processOrView == "process") 
    {
      blankTab_table <- cbind(parameters[["categoricalDataAndTime"]][ , "Sample Name"], liveReplaceBlkTable()[[1]])
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
    
    rowReplacementIndex = rowIndexInMain$custom[[input$sliderInput_BlankTab_rowsToReplace]]
    
    req(!is.null(liveReplaceBlkTable()) && !is.null(rowReplacementIndex))
    
    process$blk_ratio[[1]][rowReplacementIndex, ] <- liveReplaceBlkTable()[[1]][rowReplacementIndex, , drop = FALSE]
    process$blk_ratio[[2]][rowReplacementIndex, ] <- liveReplaceBlkTable()[[2]][rowReplacementIndex, , drop = FALSE]
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
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(rowIndexInMain$drift) | !is.integer(input$driftTab_numericInput_analyteNumber)){return()}
    driftSignal <- process$ratio_cor_b()[[1]]
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
  
  yplus <- reactive({(process$ratio_cor_b()[[1]] + process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[rowIndexInMain$drift, input$driftTab_numericInput_analyteNumber]})
  yminus <- reactive({(process$ratio_cor_b()[[1]] - process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[rowIndexInMain$drift, input$driftTab_numericInput_analyteNumber]})
  
  output$driftPlot <- renderPlot({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(rowIndexInMain$drift)){return()}
    
    elementName = parameters$analyteNames[input$driftTab_numericInput_analyteNumber]
    driftTime = parameters[["categoricalDataAndTime"]][ , "Delta Time"][rowIndexInMain$drift]
    driftValue = process$ratio_cor_b()[[1]][rowIndexInMain$drift,input$driftTab_numericInput_analyteNumber]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = elementName, pch=21, bg="lightblue")

    if (parameters$driftCorrectedElements[input$driftTab_numericInput_analyteNumber] == TRUE){
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
                                    "Internal Standard matrix (adapted for analytes)" = process$internalStandardMatrixAdaptedToAnalytes()[["signal"]],
                                    "Internal Standard ratio" = process$ratio()[["signal"]],
                                    "Blank matrix" = process$blk_ratio[["signal"]],
                                    "Blank corrected signal (CPS or ratio)" = process$ratio_cor_b()[["signal"]],
                                    "Drift matrix" = process$driftFactorMatrix()[["signal"]],
                                    "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedMatrix()[["signal"]],
                                    "Concentration" = process$concentration()[["signal"]]
    )
    
    downloadTabDataRsd <- switch(input$dataReductionStateSelection,
                                 "Counts per second" = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                 "Internal Standard matrix (adapted for analytes)" = process$internalStandardMatrixAdaptedToAnalytes()[["RSD"]],
                                 "Internal Standard ratio" = process$ratio()[["RSD"]],
                                 "Blank matrix" = process$blk_ratio[["RSD"]],
                                 "Blank corrected signal (CPS or ratio)" = process$ratio_cor_b()[["RSD"]],
                                 "Drift matrix" = process$driftFactorMatrix()[["RSD"]],
                                 "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedMatrix()[["RSD"]],
                                 "Concentration" = process$concentration()[["RSD"]]
    )
    
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
      write.csv(downloadTabData(), file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  )
}
