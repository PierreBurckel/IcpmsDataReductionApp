library(DT)
library(stringr)
library(plotly)
library(chemCal)
library(propagate)
library(MASS)
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

ICPMS_server <- function(input, output, session) {
  
  # Variable initialization ---------------------------------------------------
  
  select_input_names <- c("indexISTDchoiceWhat", "indexBlkchoiceWhat", "indexISTDchoiceIn", "indexBlkchoiceIn",
                          "selectDriftIndex", "viewConcentrationIndex")
  
  data_modified <- reactiveValues()
  data_modifiers <- reactiveValues(ISTD=list(), blank=list())
  index <- reactiveValues(custom = list())
  parsed_data <- reactiveValues()
  parsed_info <- reactiveValues()
  process_parameters <- reactiveValues(driftCorrectedElement=character(), calibration=list())
  uploaded_file <- reactiveValues(main_file = NULL, std_file = NULL, ISTD_file = NULL)
  
  blank_mod_nb <- reactiveVal(0)
  extracted <- reactiveVal(0)
  ISTD_mod_nb <- reactiveVal(0)
  set_drift <- reactiveVal(0)
  
  # Reactive expressions ------------------------------------------------------
  
  # Expression holding the current user uploaded file
  temp_data_file <- reactive({input$file})
  
  
  extraction_ready <- reactive({
    !is.null(uploaded_file$main) & !is.null(uploaded_file$std)
    })
  
  is_valid_ISTD <- reactive({
    
    extracted()
    
    return(!is.null(index$ISTD) & !is.null(parsed_data$ISTD))
  })
  
  analyte <- reactive({
    
    req(extracted())
    
    return(list(value = parsed_data$element_value[ , index$analyte, drop = FALSE],
                SD = parsed_data$element_value[ , index$analyte, drop = FALSE]))
  })
  
  ISTD <- reactive({
    
    req(extracted())
    
    return(list(value = parsed_data$element_value[ , index$ISTD, drop = FALSE],
                SD = parsed_data$element_value[ , index$ISTD, drop = FALSE]))
  })
  
  parsed_data$analyte_zero_val <- reactive({
    
    req(extracted())
    
    zero_df <- matrix(0, nrow = parsed_info$smp_nb, ncol = analyteNumber())
    colnames(zero_df) <- analyteNames()
    analyte_zero_val <- list(value = zero_df, SD = zero_df)
    
    return(analyte_zero_val)
  })
  
  analyteNames <- reactive({
    
    req(extracted())
    
    return(parsed_info$element_names[index$analyte])
  })
  
  analyteNumber <- reactive({
    
    req(extracted())
    
    return(sum(index$analyte))
  })
  
  ISTDNames <- reactive({
    
    req(extracted())
    
    return(parsed_info$element_names[index$ISTD])
  })
  
  ISTDNumber <- reactive({
    
    req(extracted())
    
    return(sum(index$ISTD))
  })
  
  activeISTDmodifier <- reactive({
    
    req(extracted())

    elements <- seq(ISTDNumber())
    method <- input$ISTDcorrectionMethod
    whatIndex <- index$custom[[input$indexISTDchoiceWhat]]
    inIndex <- index$custom[[input$indexISTDchoiceIn]]
    
    return(list(elements = elements, method = method,
                whatIndex = whatIndex, inIndex = inIndex))
  })

  activeBlankModifier <- reactive({
    
    req(extracted())
    
    elements <- seq(analyteNumber())
    method <- input$blkInterpolationMethod
    whatIndex <- index$custom[[input$indexBlkchoiceWhat]]
    inIndex <- index$custom[[input$indexBlkchoiceIn]]
    
    return(list(elements = elements, method = method,
                whatIndex = whatIndex, inIndex = inIndex))
  })
  
  data_modified$ISTD <- reactive({
    
    req(extracted())
    
    return(getModifiedData(dat = ISTD(), modifiers = data_modifiers$ISTD,
                           selectedModifiers = input$ISTDmodifiers))
  })
  
  ISTDmatrix <- reactive({
    
    req(extracted())
    
    if(is.null(parsed_data$ISTD)) {
      return(matrix(1, parsed_info$smp_nb, analyteNumber()))
    } else {
      return(createISTDMatrix(parsed_data$ISTD, data_modified$ISTD()))
    }
  })
  
  data_modified$ratio <- reactive({
    
    req(extracted())
    
    return(propagateUncertainty(a = analyte(), b = ISTDmatrix(),
                                operation="division"))
  })
  
  data_modified$blankRatio <- reactive({
    
    req(extracted())
    
    if (input$useBlankCorrection == FALSE) {
      return(parsed_data$analyte_zero_val())
    } else {
      return(getModifiedData(dat = data_modified$ratio(), modifiers = data_modifiers$blank,
                             selectedModifiers = input$blankModifiers))
    }
  })

  data_modified$blankCorrectedRatio <- reactive({
    
    req(extracted())
    
    return(propagateUncertainty(a = data_modified$ratio(), b = data_modified$blankRatio(),
                                operation="substraction"))
  })
  
  driftData <- reactive({
    
    req(index$drift, data_modified$blankCorrectedRatio)
    
    return(list(value = data_modified$blankCorrectedRatio()[["value"]][index$drift, , drop=FALSE],
                SD = data_modified$blankCorrectedRatio()[["SD"]][index$drift, , drop=FALSE]))
  })
  
  # Complex reactive expressions ----------------------------------------------
  
  # R models used for external drift correction
  data_modified$driftModels <- reactive({
    
    req(driftData())
    
    driftModels <- list()
    
    for (element in analyteNames()) {

      driftTime <- as.numeric(parsed_info$delta_t[index$drift])
      driftValues <- driftData()[["value"]][ , element]

      if (process_parameters$driftCorrectedElement[[element]] == "None") {
        driftModels[[element]] <- NULL
      }
      else if (process_parameters$driftCorrectedElement[[element]] == "Linear") {
        driftModels[[element]] <- lm(driftValues ~ driftTime)
      }
      else if (process_parameters$driftCorrectedElement[[element]] == "Quadratic") {
        driftModels[[element]] <- lm(driftValues ~ poly(driftTime, 2, raw = TRUE))
      }
    }
    
    return(driftModels)
    
  })
  
  # External drift factors for all analytes computed from the drift models
  data_modified$driftFactor <- reactive({
    
    req(index$drift)
    
    driftFactor <- list(value = vector(), SD = vector()) # drift factor for all elements
    analyteDriftFactor <- list(value = vector(), SD = vector()) # drift factor for analyte being evaluated
    
    delta_t <- as.numeric(parsed_info$delta_t)
    
    for (element in analyteNames()) {
      
      # Drift evaluation starts after the calibration of a given element
      driftStart <- getElementDriftIndex(element, parsed_data$std, parsed_info$lvl_col, index$drift)[1]
      drift_nb <- parsed_info$smp_nb - driftStart + 1
      non_drift_nb <- parsed_info$smp_nb - drift_nb
      
      analyteDriftFactor[["value"]] <- rep(1, drift_nb)
      analyteDriftFactor[["SD"]] <- rep(0, drift_nb)
      
      t0 <- delta_t[driftStart]
      
      driftModel <- data_modified$driftModels()[[element]]
      
      if (!is.null(driftModel)) {
        
        # Drift prediction and conversion to relative
        
        timeInterval <- delta_t[delta_t >= t0]
        
        driftPredict <- predict(driftModel, newdata = data.frame(driftTime = timeInterval),
                                interval="predict", confidence = 0.68)
        prediction_val <- list(value = driftPredict[ , "fit"],
                               SD = (driftPredict[ , "upr"] - driftPredict[ , "fit"]))
        first_prediction_val <- (list(value = prediction_val[["value"]][1],
                                      SD = prediction_val[["SD"]][1]))
        
        analyteDriftFactor <- propagateUncertainty(a=prediction_val, b=first_prediction_val,
                                                   operation="division")
        # First relative drift value is of 1 and SD of 0
        analyteDriftFactor[["value"]][1, ] <- 1
        analyteDriftFactor[["SD"]][1, ] <- 0
      }
      
      # Concatenation of values drift- and non-drift corrected
      # Creation of the driftFactor matrix containing relative drift values for all analytes
      
      driftFactor[["value"]] <- cbind(driftFactor[["value"]],
                                      c(rep(1, non_drift_nb), analyteDriftFactor[["value"]]))
      driftFactor[["SD"]] <- cbind(driftFactor[["SD"]],
                                   c(rep(0,driftStart-1),analyteDriftFactor[["SD"]]))
      
    }
    
    return(driftFactor)
    
  })

  # Standard signal (value and SD) of all analytes
  data_modified$standardData <- reactive({
    
    req(data_modified$blankCorrectedRatio())
    
    signal <- data_modified$blankCorrectedRatio()

    standardData <- list()
    
    for (element in analyteNames()) {
      calibration_data <- getCalibrationData(elementFullName = element, signal=signal,
                                            stdIdentificationColumn=parsed_info$lvl_col, stdDataFrame = parsed_data$std)
      standardData[[element]] <- calibration_data
    }
    
    return(standardData)
    
  })

  # Matrix format calibration models of all analytes
  data_modified$calibration_models <- reactive({
    
    req(data_modified$standardData())
    
    calibration_models <- list()
    
    for (element in analyteNames()) {
      calibration_data <- data_modified$standardData()[[element]]
      calibration_models[[element]] <- getCalibrationModel(element, process_parameters, calibration_data)
    }
    
    return(calibration_models)
    
  })

  # Analyte concentrations
  data_modified$concentration <- reactive({
    
    req(data_modified$blankCorrectedRatio())
    
    concentration_table <- vector()
    
    for (element in analyteNames()) {
      
      calibration_model <- data_modified$calibration_models()[[element]]
      blank_cor_value <- data_modified$blankCorrectedRatio()[["value"]]
      blank_cor_sd <- data_modified$blankCorrectedRatio()[["SD"]]
      drift_value <- data_modified$driftFactor()[["value"]]
      drift_sd <- data_modified$driftFactor()[["SD"]]
      
      model_parameters <- calibration_model[["modelParam"]]
      covariance_matrix <- calibration_model[["covMatrix"]]
      signal <- blank_cor_value[ , element]
      signal_variance <- blank_cor_sd[ , element]^2
      drift <- drift_value[ , element]
      drift_variance <- drift_sd[ , element]^2
      
      analyte_concentrations <- getConcentration(m = model_parameters, Cm = covariance_matrix,
                                                sig = signal, vsig = signal_variance,
                                                dft = drift, vdft = drift_variance)
      
      concentration_table <- cbind(concentration_table, analyte_concentrations)
    }
    
    return(concentration_table)
    
  })
  
  # File import ---------------------------------------------------------------
  
  # Text display of imported file
  
  output$main_assignment_txt <- renderText({
    main_file <- uploaded_file$main
    renderState(!is.null(main_file), stateTxt = "Main file name: ", invalidStateTxt = "Unassigned", validStateTxt = main_file$name)
  })
  output$std_assignment_txt <- renderText({
    std_file <- uploaded_file$std
    renderState(!is.null(std_file), stateTxt = "Standard file name: ", invalidStateTxt = "Unassigned", validStateTxt = std_file$name)
  })
  output$ISTD_assignment_txt <- renderText({
    ISTD_file <- uploaded_file$ISTD
    renderState(!is.null(ISTD_file), stateTxt = "ISTD file name: ", invalidStateTxt = "Unassigned", validStateTxt = ISTD_file$name)
  })
  
  output$extract_ready_txt <- renderText({
    renderState(extraction_ready(), stateTxt = "Data extraction ", invalidStateTxt = "impossible", validStateTxt = "ready")
  })
  
  # Table rendering of current uploaded file
  
  output$tempFilePreview <- renderTable({
    
    req(temp_data_file())
    
    temp_data_file <- read.table(temp_data_file()$datapath, sep =";", header = FALSE)
    
    rowNumber <- min(nrow(temp_data_file), input$fileUpload_nrow) 
    columnNumber <- min(length(temp_data_file), input$fileUpload_ncolumn) 
    
    temp_data_file[1:rowNumber, 1:columnNumber]
  })
  
  # Set main, std and ISTD files
  
  observeEvent(input$setAsMain, {
    req(temp_data_file())
    uploaded_file$main <- temp_data_file()
  })
  
  observeEvent(input$setAsStd, {
    req(temp_data_file())
    uploaded_file$std <- temp_data_file()
  })
  
  observeEvent(input$setAsISTD, {
    req(temp_data_file())
    uploaded_file$ISTD <- temp_data_file()
  })
  
  # Download ISTD template
  
  output$downloadISTDTemplate <- downloadHandler(
    filename = "ISTD_Template.csv",
    content = function(file) {
      write.csv(createISTDtemplate((uploaded_file$main)$datapath),
                file, sep=";", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
      })

  # Extraction ----------------------------------------------------------------
  
  observeEvent(input$extract, {
    
    req(extraction_ready())
    
    main_file <- uploaded_file$main
    std_file <- uploaded_file$std
    ISTD_file <- uploaded_file$ISTD
    
    parse_function <- get_parser("agilent")
    
    parsed <- parse_function(main_file$datapath,
                             std_file$datapath,
                             ISTD_file$datapath)
    
    parsed_data$element_value <- parsed[["element_value"]]
    parsed_data$element_sd <- parsed[["element_sd"]]
    parsed_data$std <- parsed[["std"]]
    
    parsed_info$element_names <- parsed[["element_names"]]
    parsed_info$smp_nb <- parsed[["smp_nb"]]
    parsed_info$time_col <- parsed[["time_col"]]
    parsed_info$name_col <- parsed[["name_col"]]
    parsed_info$type_col <- parsed[["type_col"]]
    parsed_info$lvl_col <- parsed[["lvl_col"]]
    parsed_info$delta_t <- parsed[["delta_t"]]
    
    index$analyte <- parsed[["analyte_index"]]
    index$ISTD <- parsed[["ISTD_index"]]
    index$numericalColumns <- parsed[["num_col_index"]]
    index$custom[["All"]] <- which(rep(x= TRUE, parsed_info$smp_nb))
    
    extracted(extracted() + 1)
    
    # If there exists an ISTD file that has been uploaded, extract and store it
    
    if (!is.null(ISTD_file)) {
      parsed_data$ISTD <- read.table(ISTD_file$datapath,
                                     header = TRUE, sep=';',
                                     stringsAsFactors=FALSE)
      
      parsed_data$ISTD <- as.matrix(parsed_data$ISTD)
    }
    
    # Initializing process parameters
    
    process_parameters$driftCorrectedElement <- rep("None", analyteNumber())
    names(process_parameters$driftCorrectedElement) <- analyteNames()
    
    process_parameters$forceIntercept <- rep(FALSE, analyteNumber())
    names(process_parameters$forceIntercept) <- analyteNames()
    
    process_parameters$autoAdaptCalibration <- rep(FALSE, analyteNumber())
    names(process_parameters$autoAdaptCalibration) <- analyteNames()

    process_parameters$useWeithedRegression  <- rep(FALSE, analyteNumber())
    names(process_parameters$useWeithedRegression) <- analyteNames()
    
    process_parameters$regressionWeight <- rep("1/var", analyteNumber())
    names(process_parameters$regressionWeight) <- analyteNames()
  })
  
  # Index creation ------------------------------------------------------------
  
  # Table
  
  output$index_table <- DT::renderDT({

    req(extracted())
    
    searchWhere = input$searchIndexwhere
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchWhere == 'smp'){
      header_col <- parsed_info$name_col
      } 
    else if (searchWhere == 'lvl'){
      header_col = parsed_info$lvl_col
      }
    else if (searchWhere == 'type'){
      header_col = parsed_info$type_col
      }
    
    if (searchType == 'ematch') searchWhat <- paste("^", searchWhat, "$", sep="")
    
    if (searchWhat == ""){
      index$temp <- seq(parsed_info$smp_nb)
    } else {
      index$temp <- grep(searchWhat, header_col)
    }
    
    if (displayWhat == "ISTD" & !is.null(ISTD()[["value"]])) {
      index_table <- ISTD()[["value"]][index$temp, , drop = FALSE]
    }
    else if (displayWhat == "analytes" & !is.null(analyte()[["value"]])) {
      index_table <- analyte()[["value"]][index$temp, , drop = FALSE]
    }
    
    row.names(index_table) <- header_col[index$temp]
    
    index_table
    
  }, options = list(dom = '', pageLength = parsed_info$smp_nb, ordering=T))
  
  # Inputs
  
  # Name of index autofilling while searching
  observeEvent(input$searchIndexwhat, {
    updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
  })
  
  # Create index button
  observeEvent(input$searchIndexCreate, {
    
    index_name <- input$searchIndexName
    custom_index <- index$temp[input$index_table_rows_selected]
    
    index$custom[[index_name]] <- custom_index
    
    updateMultipleSelectInput(session = session, input_names = select_input_names,
                              choices = names(index$custom), selected = names(index$custom)[1])
    
  }, ignoreInit = TRUE)
  
  # Table proxy creation necessary for row selection through event
  index_tableProxy <- DT::dataTableProxy("index_table", session = session)
  
  # Behavior of selectAll button
  observeEvent(input$indexSelectAll, {
    if (is.null(input$index_table_rows_selected)){
      DT::selectRows(index_tableProxy, c(1:length(index$temp)))
    }
    else if (all(input$index_table_rows_selected == seq(index$temp))){
      DT::selectRows(index_tableProxy, NULL)
    } else {
      DT::selectRows(index_tableProxy, seq(index$temp))
    }
  })
  
  # ISTD settings -------------------------------------------------------------
  
  # Table
  
  output$ISTDtable <- DT::renderDT({
    
    if (!is_valid_ISTD()) return()
    
    ISTDmode <- input$ISTDinteractionMode
    
    if (ISTDmode == "view") {
      
      df_view <- data_modified$ISTD()[["value"]]
      df_view <- format(df_view, digits = 3, scientific=T)
      
      return(df_view)
      
    }
    else if (ISTDmode == "process") {
      
      df_process <- getModifiedData(ISTD(),
                                    list(active = activeISTDmodifier()),
                                    "active")
      df_process <- ISTDprocessTable[["value"]]
      df_process <- format(df_process, digits = 3, scientific=T)
      
      return(df_process)
      
      }
  }, options = list(dom = '', pageLength = parsed_info$smp_nb, ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '150px', targets = "_all"))))
  
  # Defines a proxy for changing selections in the table
  ISTDtableProxy <- DT::dataTableProxy("ISTDtable", session = session)
  
  # Keep selection after table update
  observeEvent({
    input$indexISTDchoiceWhat
    input$indexISTDchoiceIn
    input$ISTDinteractionMode
    input$ISTDcorrectionMethod
    }, {
    if (!is.null(input$indexISTDchoiceWhat) & is_valid_ISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setISTDcorrectionMethod, {
    
    ISTD_mod_nb(ISTD_mod_nb() + 1)
    
    new_mod_id <- paste("mod", as.character(ISTD_mod_nb()), sep="")
    data_modifiers$ISTD[[new_mod_id]] <- activeISTDmodifier()
    
    mod_id <- names(data_modifiers$ISTD)
    mod_displayed_names <- sapply(data_modifiers$ISTD, getModifierName)
    names(mod_id) <- mod_displayed_names
    currentlySelected <- c(input$ISTDmodifiers, new_mod_id)
    
    updateCheckboxGroupInput(session, "ISTDmodifiers", label = "Modifiers:", choices = mod_id, selected = currentlySelected)
  })
  
  # Blank settings ------------------------------------------------------------
  
  # Disable blank process option when blank not corrected
  observeEvent(input$useBlankCorrection, {
    if (input$useBlankCorrection == FALSE) {
      updateSelectInput(session,"blkInteractionMode",
                        choices = c("View" = "view"),
                        selected= "view")
    } else {
      updateSelectInput(session,"blkInteractionMode",
                        choices = c("View" = "view", "Process" = "process"),
                        selected= "view")
    }
  }, ignoreInit = FALSE)
  
  output$blkTable <- DT::renderDT({
    
    if (is.null(analyte()[["value"]])) return()
    
    blkMode <- input$blkInteractionMode
    
    if (blkMode == "view") {
      
      df_view <- data_modified$blankRatio()[["value"]]
      df_view <- format(df_view, digits = 3, scientific=T)
      
      return(df_view)
      
    } else if (blkMode == "process") {
      
      df_process <- getModifiedData(data_modified$ratio(),
                                    list(active = activeBlankModifier()),
                                    "active")
      df_process <- df_process[["value"]]
      df_process <- format(df_process, digits = 3, scientific=T)
      
      return(df_process)
      
      }
  }, options = list(dom = '', pageLength = parsed_info$smp_nb, ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '120px', targets = "_all"))))
  
  
  #Defines a proxy for changing selections in the table
  blkTableProxy <- DT::dataTableProxy("blkTable", session = session)
  
  observeEvent({
    input$indexBlkchoiceWhat
    input$indexBlkchoiceIn
    input$blkInteractionMode
    input$blkInterpolationMethod
    }, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$setBlkInterpolationMethod, {
    
    blank_mod_nb(blank_mod_nb() + 1)
    
    new_mod_id <- paste("mod", as.character(blank_mod_nb()), sep="")
    data_modifiers$blank[[new_mod_id]] <- activeBlankModifier()
    
    mod_id <- names(data_modifiers$blank)
    mod_displayed_names <- sapply(data_modifiers$blank, getModifierName)
    names(mod_id) <- mod_displayed_names
    currentlySelected <- c(input$blankModifiers, new_mod_id)
    
    updateCheckboxGroupInput(session, "blankModifiers", label = "Modifiers:", choices = mod_id, selected = currentlySelected)
  })
  
  # Calibration settings ------------------------------------------------------
  
  observeEvent(input$setRegressionALL, {
    
    req(extracted())
    
    for (i in seq(analyteNumber())) {
      process_parameters$forceIntercept[[i]] <- input$forceIntercept
      process_parameters$autoAdaptCalibration[[i]] <- input$autoAdaptCalibration
      process_parameters$useWeithedRegression[[i]] <- input$useWeithedRegression
      process_parameters$regressionWeight[[i]] <- input$regressionWeight
    }
  }, ignoreInit = TRUE)
  
  # A serie of simple observe events to save regression parameters
  
  observeEvent(input$forceIntercept, {
    
    req(extracted())
    
    process_parameters$forceIntercept[[input$calibrationElement]] <- input$forceIntercept
  }, ignoreInit = TRUE)
  
  observeEvent(input$autoAdaptCalibration, {
    
    req(extracted())
    
    process_parameters$autoAdaptCalibration[[input$calibrationElement]] <- input$autoAdaptCalibration
  }, ignoreInit = TRUE)
  
  observeEvent(input$useWeithedRegression, {
    
    req(extracted())
    
    process_parameters$useWeithedRegression[[input$calibrationElement]] <- input$useWeithedRegression
  }, ignoreInit = TRUE)
  
  observeEvent(input$regressionWeight, {
    
    req(extracted())
    
    process_parameters$regressionWeight[[input$calibrationElement]] <- input$regressionWeight
  }, ignoreInit = TRUE)

  # UI update to display saved regression parameters associated with the selected element
  observeEvent(input$calibrationElement, {
    
    req(extracted())
    
    if (input$forceIntercept != process_parameters$forceIntercept[[input$calibrationElement]]) {
      updateCheckboxInput(session, "forceIntercept", label = "Force intercept", value = process_parameters$forceIntercept[[input$calibrationElement]])
    }
    if (input$autoAdaptCalibration != process_parameters$autoAdaptCalibration[[input$calibrationElement]]) {
      updateCheckboxInput(session, "autoAdaptCalibration", label = "Use min/max standards", value = process_parameters$autoAdaptCalibration[[input$calibrationElement]])
    }
    if (input$useWeithedRegression != process_parameters$useWeithedRegression[[input$calibrationElement]]) {
      updateCheckboxInput(session, "useWeithedRegression", label = "Use weighted linear regression", value = process_parameters$useWeithedRegression[[input$calibrationElement]])
    }
    if (input$regressionWeight != process_parameters$regressionWeight[[input$calibrationElement]]) {
      updateCheckboxInput(session, "regressionWeight", label = "Weight:", value = process_parameters$regressionWeight[[input$calibrationElement]])
    }
  }, ignoreInit = TRUE)
  
  # Calibration plot
  output$calibrationPlot <- renderPlotly({
    
    if (is.null(data_modified$standardData())) return()
    
    calibration_data <- data_modified$standardData()[[input$calibrationElement]]
    std_nb <- nrow(calibration_data)
    calibration_model <- data_modified$calibration_models()[[input$calibrationElement]]
    calibration_predict <-  getConcentration(m = calibration_model[["modelParam"]], Cm = calibration_model[["covMatrix"]],
                                            sig = calibration_data[ , "value"], vsig = calibration_data[ , "SD"]^2,
                                            dft = rep(1, std_nb), vdft = rep(0, std_nb))

    signal_residuals <- calibration_data[ ,"concentration"] - calibration_predict[ ,"value"]
    
    # To do. Add uncertainty on plot.
    # x_polygon <- c(x, rev(x))
    # y_polygon <- c(predictions[,"lwr"], rev(predictions[,"upr"]))
    
    calibrationPlot <- plot_ly()
    calibrationPlot <- add_trace(calibrationPlot, x=calibration_predict[,"value"], y=calibration_data[,"value"], type = 'scatter', mode = 'lines')
    calibrationPlot <- add_trace(calibrationPlot, x=calibration_data[,"concentration"], y=calibration_data[,"value"], type = 'scatter', mode = 'markers')

    residualPlot <- plot_ly()
    residualPlot <- add_trace(residualPlot, x=calibration_data[,"concentration"], y=signal_residuals, type = 'scatter', mode = 'markers')
    
    subplot(calibrationPlot, residualPlot, nrows = 2)

  })
  
  # Drift settings ------------------------------------------------------------
  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      set_drift(set_drift() + 1)
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- renderTable({
    req(set_drift())
    if (is.null(data_modified$blankCorrectedRatio()[[1]]) | is.null(index$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- data_modified$blankCorrectedRatio()[[1]]
    if (input$e_ind_drift < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$e_ind_drift)%%5 == 0))
    }
    uc <- min(lc + 4, analyteNumber())
    driftTable <- cbind(parsed_info$name_col[index$drift],driftSignal[index$drift,lc:uc])
    names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
    driftTable
    
  }, digits = -2)
  
  
  
  observeEvent(input$e_drift, {
    req(set_drift())
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,analyteNames(),fixed=TRUE))
  }, ignoreInit = TRUE)
  
  observeEvent(input$e_ind_drift, {
    req(set_drift())
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=analyteNames()[input$e_ind_drift])
  }, ignoreInit = TRUE)
  
  output$driftPlot <- renderPlotly({
    
    req(index$drift)
    if (is.null(data_modified$blankCorrectedRatio()[["value"]]) | is.null(index$drift)){return()}
    
    driftTime = parsed_info$delta_t[index$drift]
    driftValue = driftData()[["value"]][,input$e_ind_drift]
    driftSD = driftData()[["SD"]][,input$e_ind_drift]

    driftPlot <- plot_ly()
    driftPlot <- add_trace(driftPlot, x=driftTime, y=driftValue, type = 'scatter', mode = 'markers', error_y = list(array=driftSD, color = '#000000'))
    if (process_parameters$driftCorrectedElement[input$e_ind_drift] != "None") {
      iTime = driftTime[1]
      fTime = driftTime[length(driftTime)]
      timeInterval = seq(from=as.numeric(iTime), to=as.numeric(fTime), by=as.numeric((fTime - iTime)/100))
      driftModel <- data_modified$driftModels()[[input$e_ind_drift]]
      driftPredict=predict(driftModel, newdata = data.frame(driftTime = timeInterval),interval="predict")
      driftPlot <- add_trace(driftPlot, x=timeInterval, y=driftPredict[,"fit"], type = 'scatter', mode = 'lines')
      driftPlot <- add_annotations(driftPlot,x= 0,y= 1,xref = "paper",yref = "paper",
                                   text = paste("Adjusted R squared: ", summary(driftModel)$adj.r.squared), showarrow = F)
    }
    
    driftPlot
    #driftPlot <- add_trace(driftPlot, x=driftTime, y=driftValue, type = 'scatter', mode = 'markers', error_y = list(array=driftSD, color = '#000000'))
    # if (process_parameters$driftCorrectedElement[input$e_ind_drift] == TRUE){
    #   time_0 = driftTime[1]
    #   time_f = driftTime[length(driftTime)]
    #   timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
    #   
    #   dt = as.numeric(driftTime)
    #   driftModel <- lm(driftValue ~ poly(dt, degree=2, raw=TRUE))
    #   
    #   driftPredict=predict(driftModel, newdata = data.frame(dt = timeInterval))
    #   
    #   lines(x=timeInterval, y=driftPredict, col="red")
    # }
    
    # driftData <- subsetDfList(data_modified$blankCorrectedRatio(),index$drift, input$e_ind_drift)
    # uncertainty <- getUncertaintyInterval(driftData)
    # 
    # suppressWarnings(arrows(parsed_info$delta_t[index$drift], as.matrix(uncertainty[["lBound"]]), parsed_info$delta_t[index$drift], as.matrix(uncertainty[["uBound"]]), length=0.05, angle=90, code=3))
  })
  
  observeEvent(input$setDriftCorrection, {
    req(set_drift())
    if (is.null(process_parameters$driftCorrectedElement)){return()}
    process_parameters$driftCorrectedElement[input$e_ind_drift] <- input$driftModelSelection
  })

  # Data processing -----------------------------------------------------------
  
  output$conc <- renderTable({
    if (is.null(data_modified$concentration())){return()}
    data_modified$concentration()
    # print(1)
    # print(data_modified$concentration())
    # displayWhat = strtoi(input$viewConcentrationSwitch)
    # displayIndex = input$viewConcentrationIndex
    # 
    # if (displayWhat <= 2){displayMatrix <- data_modified$concentration()[[displayWhat]]}
    # else if (displayWhat == 3){
    #   displayMatrix <- mergeMatrixes(matrix1 = data_modified$concentration()[[1]], matrix2 = data_modified$concentration()[[2]], name1=NULL, name2="SD (%)")
    # }
    # else {}
    # print(2)
    # cbind(parsed_info$name_col, displayMatrix)[index$custom[[displayIndex]],]
  })
  
  #Button to download the list of analyte() and ISTD
  output$downloadData <- downloadHandler(
    filename = paste("data_", input$viewConcentrationIndex, ".csv", sep = ""),
    content = function(file) {
      selectedIndex = index$custom[[input$viewConcentrationIndex]]
      
      combinedConcSD <- mergeMatrixes(matrix1 = data_modified$concentration()[[1]], matrix2 = data_modified$concentration()[[2]], name1=NULL, name2="SD (%)")

      if (strtoi(input$viewConcentrationSwitch) <= 2){
        write.csv(cbind(parsed_info$name_col,data_modified$concentration()[[strtoi(input$viewConcentrationSwitch)]])[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
      else if (strtoi(input$viewConcentrationSwitch) == 3){
        write.csv(cbind(parsed_info$name_col, combinedConcSD)[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
    })
  
  
  # Observers -----------------------------------------------------------------
  
  observe({
    updateSelectInput(session,"calibrationElement",choices=analyteNames(),selected=analyteNames()[1])
    updateSelectInput(session,"e_drift",choices=analyteNames(),selected=analyteNames()[1])
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
  
  observe({
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
}
