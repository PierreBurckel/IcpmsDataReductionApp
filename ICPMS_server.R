# library(DT)
# library(stringr)
# library(plotly)
# library(shinyjs)
# source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

C_LETTER_KEYCODE <- 67

ICPMS_server <- function(input, output, session) {
  
  #Setting up lists of reactiveValues variables
  #uploadedFile contains the information about user-uploaded files
  uploadedFile <- reactiveValues()
  #extracted$data contains the list of important data after extraction
  extracted <- reactiveValues()
  #index contains all indexes serving as masks to display/process specific lines or columns of the raw data
  index <- reactiveValues()
  process <- reactiveValues()
  parameters <- reactiveValues()
  #tempDataFile stores the current file loaded by the user
  tempDataFile <- reactive({input$file})
  #extractionReady is TRUE or FALSE, depending on whether extraction can be done or not (required files assigned, variables names inputed)
  extractionReady <- reactive({!is.null(uploadedFile$raw) & !is.null(uploadedFile$std)})
  isValidISTD <- reactive({(!is.null(index$IS_CPS)) & (!is.null(extracted$data[["ISTD"]]))})
  
  #Line indexes
  #index$std <- reactive({
  #index$blk <- reactive({
  #index$drift <- reactive({
  CPS <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$CPS, drop=FALSE]})
  
  RSD <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$CPS_RSD, drop=FALSE]})
  
  IS_CPS <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$IS_CPS, drop=FALSE]})
  
  IS_CPS_RSD <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$IS_CPS_RSD, drop=FALSE]})
  
  elementNames <- reactive({
    header_1 = extracted$data[["header_1"]]
    header_1[index$CPS]})
  
  elementNumber <- reactive({
    length(which(index$CPS))})
  
  ISNames <- reactive({
    header_1 = extracted$data[["header_1"]]
    header_1[index$IS_CPS]})
  
  ISNumber <- reactive({
    length(which(index$IS_CPS))})
  
  sampleNumber <- reactive({
    nrow(extracted$data[["raw"]])})
  
  #Important columns
  timeColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    as.POSIXct(rawData[,which(header_2=="Acq. Date-Time")], format="%d/%m/%Y %H:%M")})
  
  nameColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Sample Name")]})
  
  typeColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Type")]})
  
  levelColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Level")]})
  
  dtimeColumn <- reactive({
    timeColumn() - timeColumn()[1]})
  
  ISTDmatrix <- reactive({
    if (isValidISTD()) {
      return(createISTDMatrix(extracted$data[["ISTD"]], list(signal=IS_CPS(), RSD=IS_CPS_RSD())))
    } else {
      return(list(signal = matrix(1, nrow = sampleNumber(), ncol = elementNumber()),
                  RSD = matrix(0, nrow = sampleNumber(), ncol = elementNumber())))
    }
  })
  
  process$ratio <- reactive({propagateUncertainty(a = list(signal=CPS(), RSD=RSD()), b = ISTDmatrix(), operation="division")})
  
  liveReplaceBlkTable <- reactive({
    
    replaceIndexWhat = index$custom[[input$indexBlkchoiceWhat]]
    replaceIndexIn = index$custom[[input$indexBlkchoiceIn]]
    
    replaceValues(process$ratio(), 1:elementNumber(), input$blkInterpolationMethod, replaceIndexWhat, replaceIndexIn)
    })
    
  
  process$ratio_cor_b <- reactive({
      propagateUncertainty(a=process$ratio(), b=process$blk_ratio, operation="substraction")
  })
  
  yplus <- reactive({(process$ratio_cor_b()[[1]] + process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[index$drift, input$e_ind_drift]})
  yminus <- reactive({(process$ratio_cor_b()[[1]] - process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[index$drift, input$e_ind_drift]})
  
  process$driftFactorMatrix <- reactive({
  
    for (elementIndex in 1:elementNumber()){
      
      elementFullName <- elementNames()[elementIndex]
      
      driftIndexAfterFirstStandard <- getElementDriftIndex(elementFullName = elementFullName, stdDataFrame = extracted$data[["std"]], 
                                                           stdIdentificationColumn=levelColumn(), driftIndex = index$drift)
      
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin <- timeColumn() - timeColumn()[driftIndexAfterFirstStandard][1]
      deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[deltaTimeWithFirstDriftAfterFirstStandardAsOrigin < 0] <- 0
      
      driftSignalAndDeltaTime <- cbind(Signal=process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex],
                                       dt=deltaTimeWithFirstDriftAfterFirstStandardAsOrigin[driftIndexAfterFirstStandard])  
      
      if (all(is.na(process$ratio_cor_b()[[1]][driftIndexAfterFirstStandard,elementIndex]))){
        driftPredict = rep(NA, sampleNumber())
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
  
      if (process$driftCorrectedElements[elementIndex] == TRUE){
        driftDataFrame[elementIndex] <- driftDataFrame[elementIndex] / summary(driftModel)$coefficients[1,1]
        driftDataFrame[is.na(driftDataFrame[,elementIndex]), elementIndex] <- 1
      } 
      else {
        driftDataFrame[elementIndex] <- rep(1, sampleNumber())
      }
    }
    
    return(list(signal = driftDataFrame, RSD = matrix(0, nrow = sampleNumber(), ncol = elementNumber())))
  })
  
  process$driftAndBlankCorrectedMatrix <- reactive({
    
    signalDriftCorrectedRatio <- process$ratio_cor_b()[[1]] / process$driftFactorMatrix()[["signal"]]
    signalDriftCorrectedRSD <- process$ratio_cor_b()[[2]]
    return(list(signal = signalDriftCorrectedRatio, RSD = signalDriftCorrectedRSD))
  })
  
  process$calibrationCoefficientMatrix <- reactive({
    
    for (elementIndex in 1:elementNumber()){
      
      elementFullName <- elementNames()[elementIndex]
      
      calibrationSignalUncertaintyConcentration <- getCalibrationData(elementFullName = elementFullName, signal=process$ratio_cor_b(),
                                                                      stdIdentificationColumn=levelColumn(), stdDataFrame = extracted$data[["std"]])
      
      calibrationModel <- lm(Concentration ~ 0+Signal, data=as.data.frame(calibrationSignalUncertaintyConcentration))
      
      if (elementIndex == 1){
        calibrationLinearRegressionSlope <- summary(calibrationModel)$coefficients[1,1]
      }
      else{
        calibrationLinearRegressionSlope[elementIndex] <- summary(calibrationModel)$coefficients[1,1]
      }
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
    
    selectedIndex <- index$custom[[input$viewConcentrationIndex]]
    
    downloadTabDataSignal <- switch(input$dataReductionStateSelection,
                                     "Counts per second" = CPS(),
                                     "Internal Standard matrix (adapted for analytes)" = ISTDmatrix()[["signal"]],
                                     "Internal Standard ratio" = process$ratio()[["signal"]],
                                     "Blank matrix" = process$blk_ratio[["signal"]],
                                     "Blank corrected signal (CPS or ratio)" = process$ratio_cor_b()[["signal"]],
                                     "Drift matrix" = process$driftFactorMatrix()[["signal"]],
                                     "Drift corrected signal (CPS or ratio, blank corrected)" = process$driftAndBlankCorrectedMatrix()[["signal"]],
                                     "Concentration" = process$concentration()[["signal"]]
    )
    
    downloadTabDataRsd <- switch(input$dataReductionStateSelection,
                                  "Counts per second" = RSD(),
                                  "Internal Standard matrix (adapted for analytes)" = ISTDmatrix()[["RSD"]],
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
  output$raw_assignment_txt <- renderText({
    raw_file = uploadedFile$raw
    renderState(!is.null(raw_file), stateTxt = "Raw file name: ", invalidStateTxt = "Unassigned", validStateTxt = raw_file$name)
  })
  output$std_assignment_txt <- renderText({
    std_file = uploadedFile$std
    renderState(!is.null(std_file), stateTxt = "Standard file name: ", invalidStateTxt = "Unassigned", validStateTxt = std_file$name)
  })
  output$ISTD_assignment_txt <- renderText({
    ISTD_file = uploadedFile$ISTD
    renderState(!is.null(ISTD_file), stateTxt = "ISTD file name: ", invalidStateTxt = "Unassigned", validStateTxt = ISTD_file$name)
  })
  
  #Table rendering of current uploaded file
  output$table <- renderTable({
    req(tempDataFile())
    tempDataTable = read.table(tempDataFile()$datapath, sep =";", header = FALSE)
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
    #Defines name space
    raw_file = uploadedFile$raw
    std_file = uploadedFile$std
    ISTD_file = uploadedFile$ISTD
    #Allows extraction if extraction condition are met (see extractionReady reactive value)
    if (!extractionReady()){
      return()
    } else {}
    
    extracted$data <-  extractData(raw_file$datapath, std_file$datapath)
    
    index$CPS <- extracted$data[["header_2"]] == "CPS" & !grepl("ISTD", extracted$data[["header_1"]])
    index$CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & !grepl("ISTD", extracted$data[["header_1"]])
    index$IS_CPS <- extracted$data[["header_2"]] == "CPS" & grepl("ISTD", extracted$data[["header_1"]])
    index$IS_CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & grepl("ISTD", extracted$data[["header_1"]])
    index$numericalColumns <- min(which(index$CPS)):length(extracted$data[["raw"]])
    
    extracted$data[["raw"]][,index$numericalColumns] <- sapply(extracted$data[["raw"]][,index$numericalColumns], as.numeric)
    
    ##If there exists an ISTD file that has been uploaded, extract and store it
    if (!is.null(ISTD_file)) {
      extracted$data[["ISTD"]]  <- read.table(ISTD_file$datapath, header = TRUE, sep=';', stringsAsFactors=FALSE)
    } else {}
    
    ##Defines an a priori ISTD variable that will be possible to modified in the next pane
    if (isValidISTD()) {
      process$blk_ratio <- propagateUncertainty(a = list(signal=CPS(), RSD=RSD()), b = ISTDmatrix(), operation="division")
    } else {
      process$blk_ratio <- list(signal=CPS(),RSD=RSD())
    }
    
    colnames(process$blk_ratio[["signal"]]) <- elementNames()
    colnames(process$blk_ratio[["RSD"]]) <- elementNames()
    
    index$custom[["All"]] <- which(rep(x= TRUE, sampleNumber()))
    
    ##Creates a boolean vector containing the decision of whether or not to correct signal drift for each element
    process$driftCorrectedElements <- rep(FALSE,elementNumber())
    
    parameters[["categoricalDataAndTime"]] <- cbind(timeColumn(), nameColumn(), typeColumn(), levelColumn(), dtimeColumn())
    colnames(parameters[["categoricalDataAndTime"]]) <- c("Time", "Sample Name", "Type", "Level", "Delta-Time")
    parameters[["elementNames"]] <- elementNames()
    parameters[["elementNumber"]] <- elementNumber()
    parameters[["internalStandardNames"]] <- ISNames()
    parameters[["internalStandardNumber"]] <- ISNumber()
    parameters[["sampleNumber"]] <- sampleNumber()
  })
  
  #Here we render warnings texts to help the user
  output$extract_ready_txt <- renderText({
    renderState(extractionReady(), stateTxt = "Data extraction ", invalidStateTxt = "impossible", validStateTxt = "ready")
  })
  output$ISTD_not_extracted_txt <- renderText({
    renderState(!(!is.null(uploadedFile$ISTD) & !is.null(extracted$data) & (isValidISTD() == FALSE)), stateTxt = "", invalidStateTxt = "Caution, ISTD not extracted", validStateTxt = NULL)
  })
  
# Index creation ----------------------------------------------------------

  output$indexTable <- DT::renderDT(datatable({

    if (is.null(CPS())){return()}
    
    searchWhere = input$searchIndexwhere
    firstColumnName = ""
    
    if (searchWhere == 'smp'){
      headerRows = nameColumn()
      firstColumnName = "Sample Names"}
    else if (searchWhere == 'lvl'){
      headerRows = levelColumn()
      firstColumnName = "Levels"}
    else if (searchWhere == 'type'){
      headerRows = typeColumn()
      firstColumnName = "Type"}
    
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchType == 'ematch'){searchWhat = paste("^", searchWhat, "$", sep="")}
    
    if (searchWhat == ""){
      index$temp <- c(1:sampleNumber())
    }
    else{
      index$temp <- grep(searchWhat, headerRows)
    }
    
    if (displayWhat == "ISTD" & !is.null(IS_CPS())){
      indexTable = cbind(headerRows[index$temp], IS_CPS()[index$temp,])
      colnames(indexTable) <- c(firstColumnName, ISNames())
    }
    else if (displayWhat == "analytes" & !is.null(CPS())){
      indexTable = cbind(headerRows[index$temp], CPS()[index$temp,])
      colnames(indexTable) <- c(firstColumnName, elementNames())
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
    customNumIndex = c(1:sampleNumber())[index$temp][sort(input$indexTable_rows_selected)]
    if (is.null(index$custom)) {
      index$custom <- list()
    }
    index$custom[[input$searchIndexName]] <- customNumIndex
  })
  
  indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
  
  observeEvent(input$indexSelectAll, {
    if (is.null(sort(input$indexTable_rows_selected))){
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
    else if (sort(input$indexTable_rows_selected) == c(1:length(index$temp))){
      DT::selectRows(indexTableProxy, NULL)
    }
    else{
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
  })
  
  observeEvent(index$custom, {
    if (is.null(index$custom)){return()}
    updateSelectInput(session,"indexBlkchoiceWhat", label  = "Replace what:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceIn", label  = "Replace in:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"selectDriftIndex", label  = "Define drift index:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(index$custom)),"All")
  })

# Blank verif/process -----------------------------------------------------

  #Render ISTD table if all conditions are met
  output$blkTable <- DT::renderDT(datatable({
    if (is.null(CPS())){return()}
    blkMode = input$blkInteractionMode
    
    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]
    
    df_view <- cbind(nameColumn(), process$blk_ratio[[1]])
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    
    df_process <- cbind(nameColumn(), liveReplaceBlkTable()[[1]][, lc:uc, drop = FALSE])
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
    repIndex = index$custom[[input$indexBlkchoiceWhat]]
    if (is.null(liveReplaceBlkTable())| is.null(repIndex)) {return()}

    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]

    process$blk_ratio[[1]][repIndex,lc:uc] <- liveReplaceBlkTable()[[1]][repIndex, lc:uc, drop = FALSE]
    process$blk_ratio[[2]][repIndex,lc:uc] <- liveReplaceBlkTable()[[2]][repIndex, lc:uc, drop = FALSE]
  })
  
  #Updates the slider for row and col selections when modifications are made in IS_CPS(), i.e. when the extract button is hit
  observeEvent(CPS(), {
    if (is.null(CPS())){return()}
    updateSliderInput(session,"blkRowSlider", max=sampleNumber(), value = c(1,sampleNumber()))
    updateSliderInput(session,"blkColSlider", max=elementNumber(), value = c(1,elementNumber()))
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
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- DT::renderDT(datatable({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(index$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- process$ratio_cor_b()[[1]]
    if (input$e_ind_drift < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$e_ind_drift)%%5 == 0))
    }
    uc <- min(lc + 4, elementNumber())
    driftTable <- cbind(nameColumn()[index$drift],driftSignal[index$drift,lc:uc])
    names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
    driftTable
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  
  observeEvent(elementNames(), {
    if (is.null(elementNames())){return()}
    updateSelectInput(session,"e_drift",choices=elementNames(),selected=elementNames()[1])
    updateNumericInput(session,"e_ind_drift", "Element number", 1, min = 1, max = length(elementNames()))
  })
  
  observeEvent(input$changeElementDisplayedForDrift, {
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,elementNames(),fixed=TRUE))
  })
  
  observeEvent(input$e_ind_drift, {
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=elementNames()[input$e_ind_drift])
  })
  
  output$driftPlot <- renderPlot({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(index$drift)){return()}
    elementName = elementNames()[input$e_ind_drift]
    driftTime = dtimeColumn()[index$drift]
    driftValue = process$ratio_cor_b()[[1]][index$drift,input$e_ind_drift]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = elementName, pch=21, bg="lightblue")

    if (process$driftCorrectedElements[input$e_ind_drift] == TRUE){
      time_0 = driftTime[1]
      time_f = driftTime[length(driftTime)]
      timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
      
      dt = as.numeric(driftTime)
      driftModel <- lm(driftValue ~ poly(dt, degree=2, raw=TRUE))
      
      driftPredict=predict(driftModel, newdata = data.frame(dt = timeInterval))
      
      lines(x=timeInterval, y=driftPredict, col="red")
    }
    arrows(dtimeColumn()[index$drift], yminus(), dtimeColumn()[index$drift], yplus(), length=0.05, angle=90, code=3)
  })
  
  observeEvent(input$setDriftCorrection, {
    if (is.null(process$driftCorrectedElements)){return()}
    process$driftCorrectedElements[input$e_ind_drift] <- !process$driftCorrectedElements[input$e_ind_drift]
  })
  
  observeEvent(input$pressedKey, {
    if (is.null(process$driftCorrectedElements) || input$pressedKeyId != C_LETTER_KEYCODE || input$tagId != "e_ind_drift"){return()}
    process$driftCorrectedElements[input$e_ind_drift] <- !process$driftCorrectedElements[input$e_ind_drift]
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
