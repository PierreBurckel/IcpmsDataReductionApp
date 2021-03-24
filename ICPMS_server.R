library(DT)
library(stringr)
library(plotly)
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_R6.R')



ICPMS_server <- function(input, output, session) {
  
  icpDataReactive <- IcpDataReactive$new( list(countValues = NULL, countSd = NULL) )
  #Setting up lists of reactiveValues variables
  #uploadedFile contains the information about user-uploaded files
  uploadedFiles <- reactiveValues()
  #extracted$data contains the list of important data after extraction
  extracted <- reactiveValues()
  #index contains all indexes serving as masks to display/process specific lines or columns of the raw data
  index <- reactiveValues()
  elementCorrections <- reactiveValues(data = list())
  process <- reactiveValues()
  #tempDataFile stores the current file loaded by the user
  tempDataFile <- reactive({input$file})
  #extractionReady is TRUE or FALSE, depending on whether extraction can be done or not (required files assigned, variables names inputed)
  extractionReady <- reactive({!is.null(uploadedFiles$raw) & !is.null(uploadedFiles$std)})
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
    replaceValues(propagateUncertainty(a = list(signal=CPS(), RSD=RSD()), b = ISTDmatrix(), operation="division"),
                  1:elementNumber(), input$blkInterpolationMethod, replaceIndexWhat, replaceIndexIn)})
  
  process$ratio_cor_b <- reactive({propagateUncertainty(a=process$ratio(), b=process$blk_ratio, operation="substraction")})
  
  process$ratio_cor_b_e <- reactive({
    
    browser()
    
    req(elementCorrections$data)
    
    ratio_cor_b_e <- process$ratio_cor_b()
    for (s in names(elementCorrections$data)) {
      ratio_cor_b_e <- propagateUncertainty(a=ratio_cor_b_e, b=elementCorrections$data[[s]], operation="substraction")
    }
    return(ratio_cor_b_e)
  })
  
  yplus <- reactive({(process$ratio_cor_b_e()[[1]] + process$ratio_cor_b_e()[[2]]/100*process$ratio_cor_b_e()[[1]])[index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]})
  yminus <- reactive({(process$ratio_cor_b_e()[[1]] - process$ratio_cor_b_e()[[2]]/100*process$ratio_cor_b_e()[[1]])[index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]})
  ###########################File import and viewing###########################
  
  #Text display of imported file
  output$raw_assignment_txt <- renderText({
    raw_file = uploadedFiles$raw
    renderState(!is.null(raw_file), stateTxt = "Raw file name: ", invalidStateTxt = "Unassigned", validStateTxt = raw_file$name)
  })
  output$std_assignment_txt <- renderText({
    std_file = uploadedFiles$std
    renderState(!is.null(std_file), stateTxt = "Standard file name: ", invalidStateTxt = "Unassigned", validStateTxt = std_file$name)
  })
  output$ISTD_assignment_txt <- renderText({
    ISTD_file = uploadedFiles$ISTD
    renderState(!is.null(ISTD_file), stateTxt = "ISTD file name: ", invalidStateTxt = "Unassigned", validStateTxt = ISTD_file$name)
  })
  
  #Table rendering of current uploaded file
  output$countParseVizualisationTable <- DT::renderDT(datatable({
    
    req( icpDataReactive$getCountValues() )
    req( input$countParseChoiceOfMetadata )
    
    displayedData <- switch(   
      input$countParseChoiceOfCountVizualisation,   
      "countValues"= icpDataReactive$getCountValues(),   
      "countSd"= icpDataReactive$getCountSd(),   
      "countRsd"= icpDataReactive$getCountRsd()
    )
    
    displayedDataRowNames <- icpDataReactive$getMetadata()[ , input$countParseChoiceOfMetadata]
    row.names(displayedData) <- displayedDataRowNames
    
    
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  #Buttons to set raw, std and ISTD files
  observeEvent(input$setAsRaw, {
    req(tempDataFile())
    uploadedFiles$raw <- tempDataFile()
  })
  observeEvent(input$setAsStd, {
    req(tempDataFile())
    uploadedFiles$std <- tempDataFile()
  })
  observeEvent(input$setAsISTD, {
    req(tempDataFile())
    uploadedFiles$ISTD <- tempDataFile()
  })
  
  #Button to download the list of analyte and ISTD
  output$downloadISTDTemplate <- downloadHandler(
    filename = "ISTD_Template.csv",
    content = function(file) {
      write.csv(createISTDtemplate((uploadedFiles$raw)$datapath),
                file, sep=";", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
      })

  #Button to extract important information and signal of the dataframe
  observeEvent(input$extract, {
    #Defines name space
    browser()
    rawFile = uploadedFiles$raw
    # std_file = uploadedFiles$std
    # ISTD_file = uploadedFiles$ISTD
    #Allows extraction if extraction condition are met (see extractionReady reactive value)
    countValuesAndRsd <- ExtractCountValuesAndRsdFromRawAgilentFile( rawFile$datapath )
    icpMetadata <- GetMetadataFromRawAgilentFile( rawFile$datapath )
    
    icpDataReactive$setCountValuesAndRsd( countValuesAndRsd )
    # UpdateElementNamesUiDependencies( session, icpDataReactive$getElementFullNames() )
    
    icpDataReactive$setMetadata( icpMetadata )
    updateSelectInput(session,"countParseChoiceOfMetadata", label = "Row names :",
                      choices=icpDataReactive$getMetadataNames(),selected=icpDataReactive$getMetadataNames()[1])
    updateSelectInput(session,"calibrationElement",choices=elementNames(),selected=elementNames()[1])
    updateSelectInput(session,"e_drift",choices=elementNames(),selected=elementNames()[1])
    updateSelectInput(session,"interferedElement",choices=elementNames(),selected=elementNames()[1])
    updateSelectInput(session,"interferingElement",choices=elementNames(),selected=elementNames()[1])
    # UpdateMetadataUiDependencies( session, icpDataReactive$getMetadata() )
    
    
    
    # index$CPS <- extracted$data[["header_2"]] == "CPS" & !grepl("ISTD", extracted$data[["header_1"]])
    # index$CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & !grepl("ISTD", extracted$data[["header_1"]])
    # index$IS_CPS <- extracted$data[["header_2"]] == "CPS" & grepl("ISTD", extracted$data[["header_1"]])
    # index$IS_CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & grepl("ISTD", extracted$data[["header_1"]])
    # index$numericalColumns <- min(which(index$CPS)):length(extracted$data[["raw"]])
    # 
    # extracted$data[["raw"]][,index$numericalColumns] <- sapply(extracted$data[["raw"]][,index$numericalColumns], as.numeric)
    # 
    # ##If there exists an ISTD file that has been uploaded, extract and store it
    # if (!is.null(ISTD_file)) {
    #   extracted$data[["ISTD"]]  <- read.table(ISTD_file$datapath, header = TRUE, sep=';', stringsAsFactors=FALSE)
    # } else {}
    # 
    # ##Defines an a priori ISTD variable that will be possible to modified in the next pane
    # if (isValidISTD()) {
    #   process$blk_ratio <- propagateUncertainty(a = list(signal=CPS(), RSD=RSD()), b = ISTDmatrix(), operation="division")
    # } else {
    #   process$blk_ratio <- list(signal=CPS(),RSD=RSD())
    # }
    # 
    # index$custom[["All"]] <- which(rep(x= TRUE, sampleNumber()))
    # 
    # ##Creates a boolean vector containing the decision of whether or not to correct signal drift for each element
    # process$driftCorrectedElements <- rep(FALSE,elementNumber())
    
  })
  
  #Here we render warnings texts to help the user
  output$extract_ready_txt <- renderText({
    renderState(extractionReady(), stateTxt = "Data extraction ", invalidStateTxt = "impossible", validStateTxt = "ready")
  })
  output$ISTD_not_extracted_txt <- renderText({
    renderState(!(!is.null(uploadedFiles$ISTD) & !is.null(extracted$data) & (isValidISTD() == FALSE)), stateTxt = "", invalidStateTxt = "Caution, ISTD not extracted", validStateTxt = NULL)
  })
  
  ##################################Index creation######################
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
    customNumIndex = c(1:sampleNumber())[index$temp][input$indexTable_rows_selected]
    if (is.null(index$custom)) {
      index$custom <- list()
    }
    index$custom[[input$searchIndexName]] <- customNumIndex
  })
  
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
  
  observeEvent(index$custom, {
    if (is.null(index$custom)){return()}
    updateSelectInput(session,"indexBlkchoiceWhat", label  = "Replace what:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceIn", label  = "Replace in:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"selectDriftIndex", label  = "Define drift index:", choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(index$custom)),"All")
    updateSelectInput(session,"correctionIndex", label  = "Correction Index", choices=names(index$custom),names(index$custom)[1])
  })

##################################Blank verif/process######################
  
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
  
  ################################## Element correction ############################
  
  observeEvent(input$setCorrection, {
    
    browser()
    
    correctionIndex = index$custom[[input$correctionIndex]]
    interferedElement <- input$interferedElement
    interferingElement <- input$interferingElement
    
    ratio_cor_b_values <- process$ratio_cor_b()[[1]]
    ratio_cor_b_RSD <- process$ratio_cor_b()[[2]]
    
    if (length(correctionIndex) == 1) { # one sample used for correction
      correctionFactorValue <- ratio_cor_b_values[correctionIndex, interferedElement] / ratio_cor_b_values[correctionIndex, interferingElement]
      correctionFactorRSD <- sqrt((ratio_cor_b_RSD[correctionIndex, interferedElement]/100)^2 + (ratio_cor_b_RSD[correctionIndex, interferingElement]/100)^2)*100
    }
    else if (length(correctionIndex) > 1) { # several samples used for correction
      correctionFactorValue <- mean(ratio_cor_b_values[correctionIndex, interferedElement] / ratio_cor_b_values[correctionIndex, interferingElement])
      correctionFactorRSD <- sd(ratio_cor_b_values[correctionIndex, interferedElement] / ratio_cor_b_values[correctionIndex, interferingElement]) / correctionFactorValue * 100
    }
    else {return()}
    
    zeroDf <- as.data.frame(matrix(0, nrow = sampleNumber(), ncol = elementNumber()))
    colnames(zeroDf) <- elementNames()
    
    interferenceValueDf <- zeroDf
    interferenceRSDDf <- zeroDf
    
    interferenceValue <- correctionFactorValue * ratio_cor_b_values[ , interferingElement]
    interferenceRSD <- sqrt((correctionFactorRSD/100)^2 + (ratio_cor_b_RSD[ , interferingElement]/100)^2) * 100
    
    interferenceValueDf[, interferedElement] <- interferenceValue
    interferenceRSDDf[, interferedElement] <- interferenceRSD
    
    correctionName <- paste(make.names(input$correctionIndex), make.names(interferedElement), make.names(interferingElement), sep = "_")
    
    elementCorrections$data[[correctionName]] <- list(interferenceValueDf, interferenceRSDDf)
    updateSelectInput(session,"elementCorrection", label  = "Element corrections", choices=names(elementCorrections$data), names(elementCorrections$data)[1])
  })
  
  observeEvent(input$rmCorrection, {
    
    req(input$elementCorrection)
    
    if (input$elementCorrection %in% names(elementCorrections$data)) {
      elementCorrections$data[[input$elementCorrection]] <- NULL
      updateSelectInput(session,"elementCorrection", label  = "Element corrections", choices=names(elementCorrections$data), names(elementCorrections$data)[1])
    } 
    else {
      return()
    }
    
  })
  
  output$elementCorrectionTable <- renderTable({
    
    req(elementCorrections$data)
    
    data.frame(Corrections = names(elementCorrections$data))
    
  })
  
  ################################## Calibration verification ######################
  # observeEvent(icpDataReactive$getElementFullNames(), {
  #   browser()
  #   
  #   updateSelectInput(session,"countParseChoiceOfMetadata", label = "Row names :",
  #                     choices=icpDataReactive$getElementFullNames(),selected=icpDataReactive$getElementFullNames()[1])
  #   updateSelectInput(session,"calibrationElement",choices=elementNames(),selected=elementNames()[1])
  #   updateSelectInput(session,"e_drift",choices=elementNames(),selected=elementNames()[1])
  #   updateSelectInput(session,"interferedElement",choices=elementNames(),selected=elementNames()[1])
  #   updateSelectInput(session,"interferingElement",choices=elementNames(),selected=elementNames()[1])
  # }, ignoreInit = TRUE)
  
  output$calibrationPlot <- renderPlotly({
    if (is.null(process$ratio_cor_b_e()[[1]]) | input$calibrationElement ==""){return()}
    
    eFullName <- input$calibrationElement
    
    calibrationData <- getCalibrationData(elementFullName = eFullName, signal=process$ratio_cor_b_e(),
                       stdIdentificationColumn=levelColumn(), stdDataFrame = extracted$data[["std"]])
    
    if (input$useWeithedRegression == TRUE) {
      calibrationWeights <- getWeights(calibrationData = calibrationData, fn = input$regressionWeight)
      calibrationModel <- lm(Signal ~ 0+Concentration, data=as.data.frame(calibrationData),
                             weights = calibrationWeights)
      print(calibrationData)
    } else{
      calibrationModel <- lm(Signal ~ 0+Concentration, data=as.data.frame(calibrationData))
    }
    
    
    concentrationInterval <- c(0, calibrationData[nrow(calibrationData),"Concentration"])
    
    calibrationPredict = summary(calibrationModel)$coefficients[1,1]*concentrationInterval
    
    if (input$useWeithedRegression == TRUE) {
      signalResiduals <- ((calibrationData[,"Signal"] - summary(calibrationModel)$coefficients[1,1]*calibrationData[,"Concentration"]) * calibrationWeights)
      print(calibrationData)
    } else{
      signalResiduals <- ((calibrationData[,"Signal"] - summary(calibrationModel)$coefficients[1,1]*calibrationData[,"Concentration"]))
    }
    
    calibrationPlot <- plot_ly()
    calibrationPlot <- add_trace(calibrationPlot, x=concentrationInterval, y=calibrationPredict, type = 'scatter', mode = 'lines')
    calibrationPlot <- add_trace(calibrationPlot, x=calibrationData[,"Concentration"], y=calibrationData[,"Signal"], type = 'scatter', mode = 'markers') 
    
    residualPlot <- plot_ly()
    residualPlot <- add_trace(residualPlot, x=calibrationData[,"Concentration"], y=signalResiduals, type = 'scatter', mode = 'markers')
    
    subplot(calibrationPlot, residualPlot, nrows = 2)

  })
  
  ##################################Drift verif/process######################
  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- DT::renderDT(datatable({
    if (is.null(process$ratio_cor_b_e()[[1]]) | is.null(index$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- process$ratio_cor_b_e()[[1]]
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
  
  
  
  observeEvent(input$e_drift, {
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,elementNames(),fixed=TRUE))
  })
  
  observeEvent(input$e_ind_drift, {
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=elementNames()[input$e_ind_drift])
  })
  
  output$driftPlot <- renderPlot({
    if (is.null(process$ratio_cor_b_e()[[1]]) | is.null(index$drift)){return()}
    driftTime = dtimeColumn()[index$drift]
    driftValue = process$ratio_cor_b_e()[[1]][index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = input$e_drift, pch=21, bg="lightblue")

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
  
  output$test <- renderTable({
    if (is.null(process$driftCorrectedElements)){return()}
    process$driftCorrectedElements
  })
  ####################Process
  observeEvent(input$process, {
    if (is.null(process$driftCorrectedElements)){return()}
    process$conc <- processData(process$ratio_cor_b_e(), elementNames(), extracted$data[["std"]], index$drift,levelColumn(), timeColumn(), process$driftCorrectedElements)
  })
  
  output$conc <- DT::renderDT(datatable({
    if (is.null(process$conc)){return()}
    
    displayWhat = strtoi(input$viewConcentrationSwitch)
    displayIndex = input$viewConcentrationIndex
    
    if (displayWhat <= 2){displayMatrix <- process$conc[[displayWhat]]}
    else if (displayWhat == 3){
      displayMatrix <- mergeMatrixes(matrix1 = process$conc[[1]], matrix2 = process$conc[[2]], name1=NULL, name2="RSD (%)")
    }
    else {}
    
    cbind(nameColumn(), displayMatrix)[index$custom[[displayIndex]],]
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  # Download table
  
  output$downloadData <- downloadHandler(
    filename = paste("data_", input$viewConcentrationIndex, ".csv", sep = ""),
    content = function(file) {
      
      selectedIndex <- index$custom[[input$viewConcentrationIndex]]
      displayWhat <- strtoi(input$viewConcentrationSwitch)
      
      combinedConcRSD <- mergeMatrixes(matrix1 = process$conc[[1]],
                                      matrix2 = process$conc[[2]],
                                      name1=NULL,
                                      name2="RSD (%)")
      
      if (displayWhat <= 2){
        write.csv(cbind(nameColumn(), process$conc[[displayWhat]])[selectedIndex, ],
                  file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
      else if (displayWhat == 3){
        write.csv(cbind(nameColumn(), combinedConcRSD)[selectedIndex,],
                  file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    })
}
