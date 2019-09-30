library(DT)
library(stringr)
library(plotly)
library(MASS)
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

ICPMS_server <- function(input, output, session) {
  
  dataList <- list()
  blankList <- list()

  uploadedFile <- reactiveValues()
  index <- reactiveValues()
  dataModified <- reactiveValues()
  dataModifiers <- reactiveValues(ISTD=list(), blank=list())
  processParameters <- reactiveValues(driftCorrectedElement=logical(), calibration=list())
  
  ISTDmodNumber <- reactiveVal(0)
  blankModNumber <- reactiveVal(0)
  extracted <- reactiveVal(0)
  setDrift <- reactiveVal(0)
  
  tempDataFile <- reactive({input$file})
  
  #faire une fonction isValidRaw et isValidStd qui retournent FALSE si NULL et si non valide
  extractionReady <- reactive({!is.null(uploadedFile$raw) & !is.null(uploadedFile$std)})
  #remplacer par ISTDReady et mettre une fonction isValidISTD pour check validité
  isValidISTD <- reactive({
    extracted()
    return(!is.null(index$ISTD[["value"]]) & !is.null(dataList[["ISTD"]]))
  })
  
  analyte <- reactive({
    req(extracted())
    rawData <- dataList[["raw"]]
    return(list(value=rawData[,index$analyte[["value"]], drop=FALSE],
                SD=rawData[,index$analyte[["RSD"]], drop=FALSE]/100*rawData[,index$analyte[["value"]], drop=FALSE]))
  })

  ISTD <- reactive({
    req(extracted())
    rawData <- dataList[["raw"]]
    return(list(value=rawData[,index$ISTD[["value"]], drop=FALSE],
                SD=rawData[,index$ISTD[["RSD"]], drop=FALSE]/100*rawData[,index$ISTD[["value"]], drop=FALSE]))
  })
  
  analyteNames <- reactive({
    req(extracted())
    header1 <- dataList[["header_1"]]
    return(header1[index$analyte[["value"]]])
  })
  
  analyteNumber <- reactive({
    req(extracted())
    return(length(which(index$analyte[["value"]])))
  })
  
  ISTDNames <- reactive({
    req(extracted())
    header1 <- dataList[["header_1"]]
    return(header1[index$ISTD[["value"]]])
  })
  
  ISTDNumber <- reactive({
    req(extracted())
    return(length(which(index$ISTD[["value"]])))
  })
  
  activeISTDmodifier <- reactive({
    req(extracted())
    req(isValidISTD())
    elements <- 1:ISTDNumber()
    method <- input$ISTDcorrectionMethod
    whatIndex <- index$custom[[input$indexISTDchoiceWhat]]
    inIndex <- index$custom[[input$indexISTDchoiceIn]]
    return(list(elements=elements, method=method, whatIndex=whatIndex, inIndex=inIndex))
  })

  activeBlankModifier <- reactive({
    req(extracted())
    elements <- 1:analyteNumber()
    method <- input$blkInterpolationMethod
    whatIndex <- index$custom[[input$indexBlkchoiceWhat]]
    inIndex <- index$custom[[input$indexBlkchoiceIn]]
    return(list(elements=elements, method=method, whatIndex=whatIndex, inIndex=inIndex))
  })
  
  dataModified$ISTD <- reactive({
    req(extracted())
    return(getModifiedData(ISTD(), dataModifiers$ISTD, input$ISTDmodifiers))
  })
  
  ISTDmatrix <- reactive({
    req(extracted())
    if(is.null(dataList[["ISTD"]])) {
      return(matrix(1,sampleNumber,analyteNumber()))
    }
    else{
      return(createISTDMatrix(dataList[["ISTD"]], dataModified$ISTD()[["value"]]))
    }
  })
  
  dataModified$ratio <- reactive({
    req(extracted())
    return(list(value=analyte()[["value"]]/ISTDmatrix(), SD=analyte()[["SD"]]))
  })
  
  dataModified$blankRatio <- reactive({
    req(extracted())
    if (input$useBlankCorrection == FALSE) {
      return(blankList)
    }
    else {
      return(getModifiedData(dataModified$ratio(), dataModifiers$blank, input$blankModifiers))
    }
  })

  dataModified$blankCorrectedRatio <- reactive({
    req(extracted())
    return(propagateUncertainty(a=dataModified$ratio(), b=dataModified$blankRatio(), operation="substraction"))
  })

  dataModified$standardData <- reactive({
    req(dataModified$blankCorrectedRatio())
    signal <- dataModified$blankCorrectedRatio()
    #aac <- processParameters$autoAdaptCalibration
    standardData <- list()
    for (i in seq(analyteNumber())) {
      calibrationData <- getCalibrationData(elementFullName = analyteNames()[i], signal=signal,
                                            stdIdentificationColumn=levelColumn, stdDataFrame = dataList[["std"]])
      standardData[[analyteNames()[i]]] <- calibrationData
    }
    return(standardData)
  })

  dataModified$calibrationModels <- reactive({
    req(dataModified$standardData())
    calibrationModels <- list()
    for (i in seq(analyteNumber())) {
      element <- analyteNames()[i]
      calibrationData <- dataModified$standardData()[[element]]
      calibrationModels[[element]] <- getCalibrationModel(element, processParameters, calibrationData)
    }
    return(calibrationModels)
  })

  dataModified$concentration <- reactive({
    req(extracted())
    req(dataModified$blankCorrectedRatio())
    return(processData(dataModified$blankCorrectedRatio(), analyteNames(), dataList[["std"]], index$drift,
                                             levelColumn, timeColumn, processParameters$driftCorrectedElement, processParameters$calibration))
  })
  ###########################File import and viewing###########################
  
  #Text display of imported file
  output$raw_assignment_txt <- renderText({
    raw_file <- uploadedFile$raw
    renderState(!is.null(raw_file), stateTxt = "Raw file name: ", invalidStateTxt = "Unassigned", validStateTxt = raw_file$name)
  })
  output$std_assignment_txt <- renderText({
    std_file <- uploadedFile$std
    renderState(!is.null(std_file), stateTxt = "Standard file name: ", invalidStateTxt = "Unassigned", validStateTxt = std_file$name)
  })
  output$ISTD_assignment_txt <- renderText({
    ISTD_file <- uploadedFile$ISTD
    renderState(!is.null(ISTD_file), stateTxt = "ISTD file name: ", invalidStateTxt = "Unassigned", validStateTxt = ISTD_file$name)
  })
  
  #Table rendering of current uploaded file
  output$tempFilePreview <- renderTable({
    req(tempDataFile())
    tempDataTable <- read.table(tempDataFile()$datapath, sep =";", header = FALSE)
    rowNumber <- min(nrow(tempDataTable), input$fileUpload_nrow) 
    columnNumber <- min(length(tempDataTable), input$fileUpload_ncolumn) 
    tempDataTable[1:rowNumber, 1:columnNumber]
  })
  
  #Buttons to set raw, std and ISTD files
  observeEvent(input$setAsRaw, {
    req(tempDataFile())
    uploadedFile$raw <- tempDataFile()
  })
  observeEvent(input$setAsStd, {
    req(tempDataFile())
    uploadedFile$std <- tempDataFile()
  })
  observeEvent(input$setAsISTD, {
    req(tempDataFile())
    uploadedFile$ISTD <- tempDataFile()
  })
  
  #Button to download the list of analyte and ISTD
  output$downloadISTDTemplate <- downloadHandler(
    filename = "ISTD_Template.csv",
    content = function(file) {
      write.csv(createISTDtemplate((uploadedFile$raw)$datapath),
                file, sep=";", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
      })

  #Button to extract important information and signal of the dataframe
  observeEvent(input$extract, {
    #Allows extraction if extraction condition are met (see extractionReady reactive value)
    req(extractionReady())
    
    #Defines name space
    raw_file <- uploadedFile$raw
    std_file <- uploadedFile$std
    ISTD_file <- uploadedFile$ISTD
    
    
    dataList <<-  extractData(raw_file$datapath, std_file$datapath)
    
    index$analyte <- list(value=dataList[["header_2"]] == "CPS" & !grepl("ISTD", dataList[["header_1"]]), 
                      RSD=dataList[["header_2"]] == "CPS RSD" & !grepl("ISTD", dataList[["header_1"]]))
    index$ISTD <- list(value=dataList[["header_2"]] == "CPS" & grepl("ISTD", dataList[["header_1"]]),
                         RSD=dataList[["header_2"]] == "CPS RSD" & grepl("ISTD", dataList[["header_1"]]))
    index$numericalColumns <- min(which(index$analyte[["value"]])):length(dataList[["raw"]])
    
    dataList[["raw"]][,index$numericalColumns] <<- sapply(dataList[["raw"]][,index$numericalColumns], as.numeric)
    
    extracted(extracted() + 1)
    
    rawData <- dataList[["raw"]]
    header1 <- dataList[["header_1"]]
    header2 <- dataList[["header_2"]]
    
    sampleNumber <<- nrow(dataList[["raw"]])
    timeColumn <<- as.POSIXct(rawData[,which(header2=="Acq. Date-Time")], format="%d/%m/%Y %H:%M")
    nameColumn <<- rawData[,which(header2=="Sample Name")]
    typeColumn <<- rawData[,which(header2=="Type")]
    levelColumn <<- rawData[,which(header2=="Level")]
    dtimeColumn <<- timeColumn - timeColumn[1]
    
    ##If there exists an ISTD file that has been uploaded, extract and store it
    if (!is.null(ISTD_file)) {
      dataList[["ISTD"]]  <<- read.table(ISTD_file$datapath, header = TRUE, sep=';', stringsAsFactors=FALSE)
    }
    
    #Initialisation of some variables
    blankList <<- list(value=as.data.frame(matrix(0, nrow=sampleNumber, ncol=analyteNumber())),
                            SD=as.data.frame(matrix(0, nrow=sampleNumber, ncol=analyteNumber())))
    names(blankList[["value"]]) <<- analyteNames()
    names(blankList[["SD"]]) <<- analyteNames()
    
    index$custom[["All"]] <- which(rep(x= TRUE, sampleNumber))
    
    processParameters$driftCorrectedElement <- rep(FALSE, analyteNumber())
    names(processParameters$driftCorrectedElement) <- analyteNames()
    
    processParameters$forceIntercept <- rep(FALSE, analyteNumber())
    names(processParameters$forceIntercept) <- analyteNames()
    
    processParameters$autoAdaptCalibration <- rep(FALSE, analyteNumber())
    names(processParameters$autoAdaptCalibration) <- analyteNames()

    processParameters$useWeithedRegression  <- rep(FALSE, analyteNumber())
    names(processParameters$useWeithedRegression) <- analyteNames()
    
    processParameters$regressionWeight <- rep("1/SD", analyteNumber())
    names(processParameters$regressionWeight) <- analyteNames()
  })

  #Here we render warnings texts to help the user
  output$extract_ready_txt <- renderText({
    renderState(extractionReady(), stateTxt = "Data extraction ", invalidStateTxt = "impossible", validStateTxt = "ready")
  })
  output$ISTD_not_extracted_txt <- renderText({
    renderState(!(!is.null(uploadedFile$ISTD) & !is.null(dataList) & (isValidISTD() == FALSE)), stateTxt = "", invalidStateTxt = "Caution, ISTD not extracted", validStateTxt = NULL)
  })
  
  ##################################Index creation######################
  output$indexTable <- DT::renderDT({

    if (is.null(analyte()[["value"]])){return()}
    
    searchWhere = input$searchIndexwhere
    firstColumnName = ""
    
    if (searchWhere == 'smp'){
      headerRows = nameColumn
      firstColumnName = "Sample Names"}
    else if (searchWhere == 'lvl'){
      headerRows = levelColumn
      firstColumnName = "Levels"}
    else if (searchWhere == 'type'){
      headerRows = typeColumn
      firstColumnName = "Type"}
    
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchType == 'ematch'){searchWhat = paste("^", searchWhat, "$", sep="")}
    
    if (searchWhat == ""){
      index$temp <- c(1:sampleNumber)
    }
    else{
      index$temp <- grep(searchWhat, headerRows)
    }
    
    if (displayWhat == "ISTD" & !is.null(ISTD()[["value"]])){
      indexTable <- cbind(headerRows[index$temp], ISTD()[["value"]][index$temp,])
    }
    else if (displayWhat == "analytes" & !is.null(analyte()[["value"]])){
      indexTable <- cbind(headerRows[index$temp], analyte()[["value"]][index$temp,])
    }
    else {return()}
    
    colnames(indexTable) <- c(firstColumnName, colnames(ISTD()[["value"]])[1:length(ISTD()[["value"]])])
    
    indexTable
    
  }, options = list(dom = '', pageLength = sampleNumber, ordering=T))
  
  observeEvent(input$searchIndexwhat, {
    updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
  })
  
  observeEvent(input$searchIndexCreate, {
    indexName = input$searchIndexName
    customNumIndex = c(1:sampleNumber)[index$temp][input$indexTable_rows_selected]
    if (is.null(index$custom)) {
      index$custom <- list()
    }
    index$custom[[input$searchIndexName]] <- customNumIndex
    
    updateSelectInput(session,"indexISTDchoiceWhat",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceWhat",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexISTDchoiceIn",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceIn",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"selectDriftIndex",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex",choices=c("All", names(index$custom)),"All")
  }, ignoreInit = TRUE)
  
  indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
  
  observeEvent(input$indexSelectAll, {
    if (is.null(input$indexTable_rows_selected)){
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
    else if (input$indexTable_rows_selected == c(1:length(index$temp))){
      DT::selectRows(indexTableProxy, NULL)
    }
    else{
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
  })
  
  ##################################ISTD verif/process######################
  
  #Render ISTD table if all conditions are met
  output$ISTDtable <- DT::renderDT({
    if (!isValidISTD()){
      print("No available ISTD")
      return()
    }
    
    ISTDmode <- input$ISTDinteractionMode
    ISTDviewTable <- dataModified$ISTD()
    
    allModifiers <- dataModifiers$ISTD
    allModifiers[["active"]] <- activeISTDmodifier()
    ISTDprocessTable <- getModifiedData(ISTD(), allModifiers, c(input$ISTDmodifiers, "active"))
    
    lc <- input$ISTDcolSlider[1]
    uc <- input$ISTDcolSlider[2]
    
    df_view <- cbind(nameColumn, ISTDviewTable[["value"]])
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    
    df_process <- cbind(nameColumn, ISTDprocessTable[["value"]])
    names(df_process) <- c("Sample Name", names(df_process)[2:length(df_process)])
    df_process <- format(df_process, digits = 3, scientific=T)
    
    if (ISTDmode == "view") {df_view}
    else if (ISTDmode == "process") {df_process}
  }, options = list(dom = '', pageLength = sampleNumber, ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '150px', targets = "_all"))))
  
  #Defines a proxy for changing selections in the table
  ISTDtableProxy <- DT::dataTableProxy("ISTDtable", session = session)
  
  observeEvent(input$ISTDcolSlider, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexISTDchoiceWhat, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexISTDchoiceIn, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$ISTDinteractionMode, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$ISTDcorrectionMethod, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setISTDcorrectionMethod, {
    
    ISTDmodNumber(ISTDmodNumber() + 1)
    
    #quand on voudra pouvoir selectionner des éléments particuliers
    #lc <- input$ISTDcolSlider[1]
    #uc <- input$ISTDcolSlider[2]
    
    newModID <- paste("mod", as.character(ISTDmodNumber()), sep="")
    dataModifiers$ISTD[[newModID]] <- activeISTDmodifier()
    
    modID <- names(dataModifiers$ISTD)
    modNames <- sapply(dataModifiers$ISTD, getModifierName)
    names(modID) <- modNames
    currentlySelected <- c(input$ISTDmodifiers, newModID)
    
    #selected takes a character vector of IDs
    updateCheckboxGroupInput(session, "ISTDmodifiers", label = "Modifiers:", choices = modID, selected = currentlySelected)
  })
  
  observe({
    updateSelectInput(session,"calibrationElement",choices=analyteNames(),selected=analyteNames()[1])
    updateSelectInput(session,"e_drift",choices=analyteNames(),selected=analyteNames()[1])
    if(input$ISTDinteractionMode == "process") {
      shinyjs::enable("setISTDcorrectionMethod")
      shinyjs::enable("ISTDcolSlider")
      shinyjs::enable("ISTDcorrectionMethod")
      shinyjs::enable("indexISTDchoiceIn")
    }
    else if(input$ISTDinteractionMode == "view") {
      shinyjs::disable("setISTDcorrectionMethod")
      shinyjs::disable("ISTDcolSlider")
      shinyjs::disable("ISTDcorrectionMethod")
      shinyjs::disable("indexISTDchoiceIn")
    }
  })
  
##################################Blank verif/process######################
  #browser()
  #Render ISTD table if all conditions are met
  output$blkTable <- DT::renderDT({
    if (is.null(analyte()[["value"]])){return()}
    
    blkMode <- input$blkInteractionMode
    blankViewTable <- dataModified$blankRatio()
    
    allModifiers <- dataModifiers$blank
    allModifiers[["active"]] <- activeBlankModifier()
    blankProcessTable <- getModifiedData(dataModified$ratio(), allModifiers, c(input$blankModifiers, "active"))
    
    lc <- input$blankColSlider[1]
    uc <- input$blankColSlider[2]
    
    df_view <- cbind(nameColumn, blankViewTable[["value"]])
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    
    df_process <- cbind(nameColumn, blankProcessTable[["value"]])
    names(df_process) <- c("Sample Name", names(df_process)[2:length(df_process)])
    df_process <- format(df_process, digits = 3, scientific=T)

    if (blkMode == "view") {df_view}
    else if (blkMode == "process") {df_process}
  }, options = list(dom = '', pageLength = sampleNumber, ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '120px', targets = "_all"))))
  #browser()
  #Defines a proxy for changing selections in the table
  blkTableProxy <- DT::dataTableProxy("blkTable", session = session)
  
  
  observeEvent(input$blkColSlider, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceWhat, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceIn, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInteractionMode, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInterpolationMethod, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setBlkInterpolationMethod, {
    
    blankModNumber(blankModNumber() + 1)
    
    #quand on voudra pouvoir selectionner des éléments particuliers
    #lc <- input$ISTDcolSlider[1]
    #uc <- input$ISTDcolSlider[2]
    
    newModID <- paste("mod", as.character(blankModNumber()), sep="")
    dataModifiers$blank[[newModID]] <- activeBlankModifier()
    
    modID <- names(dataModifiers$blank)
    modNames <- sapply(dataModifiers$blank, getModifierName)
    names(modID) <- modNames
    currentlySelected <- c(input$blankModifiers, newModID)
    
    #selected takes a character vector of IDs
    updateCheckboxGroupInput(session, "blankModifiers", label = "Modifiers:", choices = modID, selected = currentlySelected)
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
  
  ################################## Calibration verification ######################
  observeEvent(input$setRegressionALL, {
    req(extracted())
    for (i in 1:analyteNumber()) {
      processParameters$forceIntercept[[i]] <- input$forceIntercept
      processParameters$autoAdaptCalibration[[i]] <- input$autoAdaptCalibration
      processParameters$useWeithedRegression[[i]] <- input$useWeithedRegression
      processParameters$regressionWeight[[i]] <- input$regressionWeight
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$forceIntercept, {
    req(extracted())
    processParameters$forceIntercept[[input$calibrationElement]] <- input$forceIntercept
  }, ignoreInit = TRUE)
  
  observeEvent(input$autoAdaptCalibration, {
    req(extracted())
    processParameters$autoAdaptCalibration[[input$calibrationElement]] <- input$autoAdaptCalibration
  }, ignoreInit = TRUE)
  
  observeEvent(input$useWeithedRegression, {
    req(extracted())
    processParameters$useWeithedRegression[[input$calibrationElement]] <- input$useWeithedRegression
  }, ignoreInit = TRUE)
  
  observeEvent(input$regressionWeight, {
    req(extracted())
    processParameters$regressionWeight[[input$calibrationElement]] <- input$regressionWeight
  }, ignoreInit = TRUE)

  observeEvent(input$calibrationElement, {
    req(extracted())
    if (input$forceIntercept != processParameters$forceIntercept[[input$calibrationElement]]) {
      updateCheckboxInput(session, "forceIntercept", label = "Force intercept", value = processParameters$forceIntercept[[input$calibrationElement]])
    }
    if (input$autoAdaptCalibration != processParameters$autoAdaptCalibration[[input$calibrationElement]]) {
      updateCheckboxInput(session, "autoAdaptCalibration", label = "Use min/max standards", value = processParameters$autoAdaptCalibration[[input$calibrationElement]])
    }
    if (input$useWeithedRegression != processParameters$useWeithedRegression[[input$calibrationElement]]) {
      updateCheckboxInput(session, "useWeithedRegression", label = "Use weighted linear regression", value = processParameters$useWeithedRegression[[input$calibrationElement]])
    }
    if (input$regressionWeight != processParameters$regressionWeight[[input$calibrationElement]]) {
      updateCheckboxInput(session, "regressionWeight", label = "Weight:", value = processParameters$regressionWeight[[input$calibrationElement]])
    }
  }, ignoreInit = TRUE)
  
  output$calibrationPlot <- renderPlotly({
    if (is.null(dataModified$standardData())){return()}
    
    calibrationData <- dataModified$standardData()[[input$calibrationElement]]
    
    calibrationModel <- dataModified$calibrationModels()[[input$calibrationElement]]
    
    concentrationInterval <- c(0, calibrationData[nrow(calibrationData),"concentration"])
    
    if (any(grepl("(Intercept)",rownames(summary(calibrationModel)$coefficients)))) {
      intercept <- summary(calibrationModel)$coefficients["(Intercept)",1]
      calibrationPredict = summary(calibrationModel)$coefficients["concentration",1]*concentrationInterval + intercept
    }
    else {
      calibrationPredict = summary(calibrationModel)$coefficients["concentration",1]*concentrationInterval
    }
    
    
    signalResiduals <- studres(calibrationModel)
    
    calibrationPlot <- plot_ly()
    calibrationPlot <- add_trace(calibrationPlot, x=concentrationInterval, y=calibrationPredict, type = 'scatter', mode = 'lines')
    calibrationPlot <- add_trace(calibrationPlot, x=calibrationData[,"concentration"], y=calibrationData[,"value"], type = 'scatter', mode = 'markers')
    
    residualPlot <- plot_ly()
    residualPlot <- add_trace(residualPlot, x=calibrationData[,"concentration"], y=signalResiduals, type = 'scatter', mode = 'markers')
    
    subplot(calibrationPlot, residualPlot, nrows = 2)

  })
  
  ##################################Drift verif/process######################
  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      setDrift(setDrift() + 1)
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- renderTable({
    req(setDrift())
    if (is.null(dataModified$blankCorrectedRatio()[[1]]) | is.null(index$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- dataModified$blankCorrectedRatio()[[1]]
    if (input$e_ind_drift < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$e_ind_drift)%%5 == 0))
    }
    uc <- min(lc + 4, analyteNumber())
    driftTable <- cbind(nameColumn[index$drift],driftSignal[index$drift,lc:uc])
    names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
    driftTable
    
  }, digits = -2)
  
  
  
  observeEvent(input$e_drift, {
    req(setDrift())
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,analyteNames(),fixed=TRUE))
  }, ignoreInit = TRUE)
  
  observeEvent(input$e_ind_drift, {
    req(setDrift())
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=analyteNames()[input$e_ind_drift])
  }, ignoreInit = TRUE)
  
  output$driftPlot <- renderPlot({
    req(setDrift())
    req(index$drift)
    if (is.null(dataModified$blankCorrectedRatio()[[1]]) | is.null(index$drift)){return()}
    driftTime = dtimeColumn[index$drift]
    driftValue = dataModified$blankCorrectedRatio()[[1]][index$drift,input$e_ind_drift]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = input$e_drift, pch=21, bg="lightblue")

    if (processParameters$driftCorrectedElement[input$e_ind_drift] == TRUE){
      time_0 = driftTime[1]
      time_f = driftTime[length(driftTime)]
      timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
      
      dt = as.numeric(driftTime)
      driftModel <- lm(driftValue ~ poly(dt, degree=2, raw=TRUE))
      
      driftPredict=predict(driftModel, newdata = data.frame(dt = timeInterval))
      
      lines(x=timeInterval, y=driftPredict, col="red")
    }
    
    driftData <- subsetDfList(dataModified$blankCorrectedRatio(),index$drift, input$e_ind_drift)
    uncertainty <- getUncertaintyInterval(driftData)
    
    suppressWarnings(arrows(dtimeColumn[index$drift], as.matrix(uncertainty[["lBound"]]), dtimeColumn[index$drift], as.matrix(uncertainty[["uBound"]]), length=0.05, angle=90, code=3))
  })
  
  observeEvent(input$setDriftCorrection, {
    req(setDrift())
    if (is.null(processParameters$driftCorrectedElement)){return()}
    processParameters$driftCorrectedElement[input$e_ind_drift] <- !processParameters$driftCorrectedElement[input$e_ind_drift]
  })

  ####################Process
  
  output$conc <- renderTable({
    if (is.null(dataModified$concentration())){return()}
    print(1)
    print(dataModified$concentration())
    displayWhat = strtoi(input$viewConcentrationSwitch)
    displayIndex = input$viewConcentrationIndex
    
    if (displayWhat <= 2){displayMatrix <- dataModified$concentration()[[displayWhat]]}
    else if (displayWhat == 3){
      displayMatrix <- mergeMatrixes(matrix1 = dataModified$concentration()[[1]], matrix2 = dataModified$concentration()[[2]], name1=NULL, name2="SD (%)")
    }
    else {}
    print(2)
    cbind(nameColumn, displayMatrix)[index$custom[[displayIndex]],]
  })
  
  #Button to download the list of analyte() and ISTD
  output$downloadData <- downloadHandler(
    filename = paste("data_", input$viewConcentrationIndex, ".csv", sep = ""),
    content = function(file) {
      selectedIndex = index$custom[[input$viewConcentrationIndex]]
      
      combinedConcSD <- mergeMatrixes(matrix1 = dataModified$concentration()[[1]], matrix2 = dataModified$concentration()[[2]], name1=NULL, name2="SD (%)")

      if (strtoi(input$viewConcentrationSwitch) <= 2){
        write.csv(cbind(nameColumn,dataModified$concentration()[[strtoi(input$viewConcentrationSwitch)]])[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
      else if (strtoi(input$viewConcentrationSwitch) == 3){
        write.csv(cbind(nameColumn, combinedConcSD)[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
    })
  
  #callModule(csvFile, "default", nameColumn = reactive(nameColumn), indexCustom = reactive(index$custom), isValidISTD = reactive(isValidISTD()), ISTD = reactive(ISTD), signal=reactive(list(ISTD,ISTD_SD)))
}
