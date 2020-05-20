library(DT)
library(shinyjs)
library(stringr)
library(plotly)
library(chemCal)
library(propagate)
library(MASS)
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')



ICPMS_server <- function(input, output, session) {

  # Variable initialization ---------------------------------------------------
  
  intervalPlot_init <- reactiveVal(0)
  
  update_interval_methods <- reactiveVal(TRUE)
  
  select_input_names <- c("indexISTDchoiceWhat", "indexBlkchoiceWhat", "indexISTDchoiceIn", "indexBlkchoiceIn",
                          "selectDriftIndex", "viewConcentrationIndex", "blkPlotIndex", "interval_index")
  expected_interval_select_input <- list(method = character(), opt1 = character(), opt2 = character())
  
  data_modified <- reactiveValues()
  data_modifiers <- reactiveValues(ISTD = list(), blank = list())
  data_tables <- reactiveValues(ISTD = matrix(), blank = matrix())
  index <- reactiveValues(custom = list())
  
  interval_tmp_mdl <- reactiveValues(id = numeric(), prev_ids = numeric(), interval_sel_inputs = list(), valid_sel_inputs = boolean(),
                                         sub_index = list(), t0 = NULL, tf = NULL, times = list(),
                                        modeled_data = list())
  interval_models <- reactiveValues()
  
  parsed_data <- reactiveValues()
  parsed_info <- reactiveValues()
  process_parameters <- reactiveValues(driftCorrectedElement = character(), calibration = list())
  uploaded_file <- reactiveValues(main_file = NULL, std_file = NULL, ISTD_file = NULL)
  
  blank_mod_nb <- reactiveVal(0)
  extracted <- reactiveVal(0)
  ISTD_mod_nb <- reactiveVal(0)
  set_drift <- reactiveVal(0)
  
  # Reactive expressions ------------------------------------------------------
  
  # Expression holding the current user uploaded file
  temp_data_file <- reactive({input$file})
  
  
  extraction_ready <- reactive({
    !is.null(uploaded_file$main) && !is.null(uploaded_file$std)
    })
  
  is_valid_ISTD <- reactive({
    
    extracted()
    
    return(!is.null(index$ISTD) && !is.null(parsed_data$ISTD))
  })
  
  # interval_select_inputs <- reactive({
  # 
  #   req(extracted())
  #   req(!is.null(input$interval_index))
  # 
  #   browser()
  # 
  #   interval_select_inputs <- list()
  #   
  #   for (id in interval_tmp_mdl$id) {
  #     
  #     interval_select_inputs[["method"]] <- c(interval_select_inputs[["method"]], input[[paste0("sel_", id)]])
  #     interval_select_inputs[["opt1"]] <- c(interval_select_inputs[["opt1"]], input[[paste0("opt_1_", id)]])
  #     interval_select_inputs[["opt2"]] <- c(interval_select_inputs[["opt2"]], input[[paste0("opt_2_", id)]])
  #     
  #   }
  #   
  #   environment(update_interval_select_inputs) <- environment()
  #   update_interval_select_inputs()
  #   
  #   return(interval_select_inputs)
  # 
  # })
  
  
    
  # valid_interval_select_inputs <- reactive({
  #   
  #   req(extracted())
  #   req(!is.null(input$interval_index))
  #   
  #   browser()
  #   
  #   if (length(interval_tmp_mdl$id) == 0) {
  #     
  #     return(FALSE)
  #     
  #   } else {
  #     
  #     for (i in seq(length(interval_tmp_mdl$id))) {
  #       
  #       id <- interval_tmp_mdl$id[i]
  #       
  #       if (is.null(input[[paste0("sel_", id)]]) ||
  #           is.null(input[[paste0("opt_1_", id)]]) ||
  #           is.null(input[[paste0("opt_2_", id)]]) ||
  #           expected_interval_select_input[["method"]][i] != input[[paste0("sel_", id)]] ||
  #           expected_interval_select_input[["opt1"]][i] != input[[paste0("opt_1_", id)]] ||
  #           expected_interval_select_input[["opt2"]][i] != input[[paste0("opt_2_", id)]]) {
  #         return(FALSE)
  #       }
  #       
  #     }
  #     
  #   }
  #   
  #   return(TRUE)
  # })
  
  interval_tmp_mdl$sub_index <- reactive({

    req(extracted())

    sub_index <- list()

    time_col <- parsed_info$time_col

    interval_index <- index$custom[[input$interval_index]]

    interval_nb <- length(interval_tmp_mdl$id)

    for (i in seq(interval_nb)) {

      t0 <- interval_tmp_mdl$slices_t0()[i]
      tf <- interval_tmp_mdl$slices_tf()[i]

      sub_index[[i]] <- time_col[interval_index] >= t0 & time_col[interval_index] <= tf

    }

    return(sub_index)

  })
  
  interval_tmp_mdl$slices_t0 <- reactive({
    
    req(extracted())
    
    slices_t0 <- NULL
    
    time_col <- parsed_info$time_col
    
    interval_index <- index$custom[[input$interval_index]]
    
    for (id in interval_tmp_mdl$id) {
      slices_t0 <- c(slices_t0, time_col[interval_index][id])
    }
    
    return(slices_t0)
    
  })
  
  interval_tmp_mdl$slices_tf <- reactive({
    
    req(extracted())
    
    slices_tf <- NULL
    
    time_col <- parsed_info$time_col
    
    sequence_tf <- tail(time_col, n = 1)
    
    interval_nb <- length(interval_tmp_mdl$id)
    
    for (i in seq(interval_nb)) {
      
      if (i < interval_nb) {
        
        slices_tf <- c(slices_tf, interval_tmp_mdl$slices_t0()[i + 1])
        
      } else if (i == interval_nb) {
        
        slices_tf <- c(slices_tf, sequence_tf)
        
      }
    }
    
    return(slices_tf)
    
  })
  
  # interval_tmp_mdl$slices_sample_times <- reactive({
  #   
  #   req(extracted())
  #   
  #   times <- list(absolute = list(), relative = list())
  #   
  #   interval_nb <- length(interval_tmp_mdl$id)
  # 
  #   for (i in seq(interval_nb)) {
  #     
  #     t0 <- interval_tmp_mdl$slices_t0()[i]
  #     tf <- interval_tmp_mdl$slices_tf()[i]
  #     
  #     index <- t0 <= parsed_info$time_col & tf >= parsed_info$time_col
  #     times[[i]] <- parsed_info$time_col[index]
  #   }
  #   
  #   return(times)
  #   
  # })
  
  # interval_tmp_mdl$modeled_data <- reactive({
  #   
  #   browser()
  #   
  #   modeled_data <- list()
  #   
  #   interval_index <- index$custom[[input$interval_index]]
  #   element <- input$intervalElement
  #   interval_nb <- length(interval_tmp_mdl$id)
  #   delta_t <- parsed_info$delta_t
  #   
  #   x <- as.numeric(delta_t[interval_index])
  #   y <- list(value = analyte()[["value"]][interval_index , element],
  #             SD = analyte()[["SD"]][interval_index , element])
  #   
  #   for (i in seq(interval_nb)) {
  #     
  #     id <- interval_tmp_mdl$id[[i]]
  #     
  #     sub_index <- interval_tmp_mdl$sub_index()[[i]]
  #   
  #     sample_time <- x[sub_index]
  #     y_sub <- list(value = y[["value"]][sub_index],
  #                   SD = y[["SD"]][sub_index])
  #     
  #     method <- input[[paste0("sel_", id)]]
  #     opt1 <- input[[paste0("opt_1_", id)]]
  #     opt2 <- input[[paste0("opt_2_", id)]]
  #     
  #     x_pred <- interval_tmp_mdl$slices_sample_times()[[i]] - parsed_info$time_col[1]
  #     measured_data_index <- as.numeric(parsed_info$delta_t) %in% x_pred
  #     y_pred <- list(value = analyte()[["value"]][measured_data_index , element],
  #                    SD = analyte()[["SD"]][measured_data_index , element])
  #     
  #     modeled_data[[i]] <- get_modeled_data(sample_time, y_sub,
  #                                           method, opt1, opt2,
  #                                           x_pred)
  #   }
  #   
  #   return(modeled_data)
  # })
  # 
  # interval_tmp_mdl$sequence_pred_val <- reactive({
  #   
  #   interval_nb <- length(interval_tmp_mdl$id)
  #   sequence_pred_val <- NULL
  #   
  #   for (i in seq(interval_nb)) {
  #     slice_pred_val <- interval_tmp_mdl$modeled_data()[[i]]
  #     if (i != interval_nb) {
  #       sequence_pred_val <- c(sequence_pred_val, slice_pred_val[-length(slice_pred_val)])
  #     } else {
  #       sequence_pred_val <- c(sequence_pred_val, slice_pred_val)
  #     }
  #   }
  #   
  #   first_slice_t0 <- interval_tmp_mdl$slices_t0()[1]
  #   non_modeled_smp_nb <- sum(parsed_info$time_col < first_slice_t0)
  #   
  #   if (non_modeled_smp_nb > 0) {
  #     sequence_pred_val <- c(sequence_pred_val, rep(sequence_pred_val[1], non_modeled_smp_nb))
  #   }
  #   
  #   return(sequence_pred_val)
  # })
  
  analyte <- reactive({
    
    req(extracted())
    
    return(list(value = parsed_data$element_value[ , index$analyte, drop = FALSE],
                SD = parsed_data$element_sd[ , index$analyte, drop = FALSE]))
  })
  
  ISTD <- reactive({
    
    req(extracted())
    
    return(list(value = parsed_data$element_value[ , index$ISTD, drop = FALSE],
                SD = parsed_data$element_sd[ , index$ISTD, drop = FALSE]))
  })
  
  parsed_data$analyte_zero_val <- reactive({
    
    req(extracted())
    
    zero_df <- matrix(0, nrow = parsed_info$smp_nb, ncol = analyteNumber())
    colnames(zero_df) <- analyteNames()
    row.names(zero_df) <- parsed_info$name_col
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
  
  # Displayed ISTD table
  data_tables$ISTD <- reactive({
    
    if (!is_valid_ISTD()) return()
    
    ISTDmode <- input$ISTDinteractionMode
    
    if (ISTDmode == "view") {
      
      df_view <- data_modified$ISTD()[["value"]]
      
      return(df_view)
      
    }
    else if (ISTDmode == "process") {
      
      df_process <- getModifiedData(ISTD(),
                                    list(active = activeISTDmodifier()),
                                    "active")
      df_process <- df_process[["value"]]
      
      return(df_process)
      
    }
  })
  
  # Displayed blank table
  data_tables$blank <- reactive({
    if (is.null(analyte()[["value"]])) return()
    
    blkMode <- input$blkInteractionMode
    blkDisplay <- input$blkDisplayMode
    
    if (blkMode == "view") {
      
      if (blkDisplay == "blank") {
        df_view <- data_modified$blankRatio()[["value"]]
      }
      else if (blkDisplay == "sig_minus_blk") {
        df_view <- data_modified$blankCorrectedRatio()[["value"]]
      }
      
      return(df_view)
      
    } else if (blkMode == "process") {
      
      df_blank <- getModifiedData(data_modified$ratio(),
                                  list(active = activeBlankModifier()),
                                  "active")[["value"]]
      
      if (blkDisplay == "blank") {
        df_process <- df_blank
      }
      else if (blkDisplay == "sig_minus_blk") {
        df_sig_minus_blk <- data_modified$ratio()[["value"]] - df_blank
        df_process <- df_sig_minus_blk
      }
      
      return(df_process)
      
    }
  })
  
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
                                interval = "predict", confidence = 0.68)
        prediction_val <- list(value = driftPredict[ , "fit"],
                               SD = (driftPredict[ , "upr"] - driftPredict[ , "fit"]))
        first_prediction_val <- (list(value = prediction_val[["value"]][1],
                                      SD = prediction_val[["SD"]][1]))
        
        analyteDriftFactor <- propagateUncertainty(a=prediction_val, b=first_prediction_val,
                                                   operation="division")
        # First relative drift value is of 1 and SD of 0
        analyteDriftFactor[["value"]][1] <- 1
        analyteDriftFactor[["SD"]][1] <- 0
      }
      
      # Concatenation of values drift- and non-drift corrected
      # Creation of the driftFactor matrix containing relative drift values for all analytes
      
      driftFactor[["value"]] <- cbind(driftFactor[["value"]],
                                      c(rep(1, non_drift_nb), analyteDriftFactor[["value"]]))
      driftFactor[["SD"]] <- cbind(driftFactor[["SD"]],
                                   c(rep(0,driftStart-1),analyteDriftFactor[["SD"]]))
      
      colnames(driftFactor[["value"]]) <- analyteNames()
      colnames(driftFactor[["SD"]]) <- analyteNames()
      
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
                                            stdIdentificationColumn=parsed_info$lvl_col,
                                            stdDataFrame = parsed_data$std)
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
  # for DT extensions see: https://rstudio.github.io/DT/extensions.html
  
  output$index_table <- DT::renderDT(datatable({

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
    
    if (displayWhat == "ISTD" && !is.null(ISTD()[["value"]])) {
      index_table <- ISTD()[["value"]][index$temp, , drop = FALSE]
    }
    else if (displayWhat == "analytes" && !is.null(analyte()[["value"]])) {
      index_table <- analyte()[["value"]][index$temp, , drop = FALSE]
    }
    
    row.names(index_table) <- header_col[index$temp]
    
    index_table
    
  }, extensions = c('Scroller', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
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
    
    # updateMultipleSelectInput(session = session, input_names = select_input_names,
    #                           choices = names(index$custom), selected = names(index$custom)[1])
    
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
  
  # Interval creation ---------------------------------------------------------
  
  observeEvent(input$create_interval, {
    
    # browser()
    
    interval_name <- input$interval_name
    
    if (interval_name == "") {
      showModal(modalDialog(
        title = "Missing name",
        "A name is required for the interval"
      ))
      return()
    }
    
    time_continuous <- NULL
    data_continuous <- NULL
    
    for (i in seq(length(interval_tmp_mdl$id))) {
      time_continuous <- c(time_continuous, interval_tmp_mdl$slices_sample_times()[[i]])
      data_continuous <- list(value = c(data_continuous[["value"]], interval_tmp_mdl$modeled_data()[["value"]][[i]]),
                              SD = c(data_continuous[["SD"]], interval_tmp_mdl$modeled_data()[["SD"]][[i]]))
    }
    
    duplicate_index <- duplicated(time_continuous)
    
    if (any(duplicate_index)) {
      time_continuous <- time_continuous[-which(duplicate_index)]
      data_continuous[["value"]] <- data_continuous[["value"]][-which(duplicate_index)]
      data_continuous[["SD"]] <- data_continuous[["SD"]][-which(duplicate_index)]
    }
    
    interval_models[[interval_name]] <- list(id = interval_tmp_mdl$id,
                                             interval_select_inputs = interval_select_inputs(),
                                             interval_index = input$interval_index,
                                             t0 = interval_tmp_mdl$slices_t0(),
                                             tf = interval_tmp_mdl$slices_tf(),
                                             times = parsed_info$time_col,
                                             modeled_data = interval_tmp_mdl$sequence_pred_val()
                                             )
    
    updateSelectInput(session, "interval_names",
                      choices = names(interval_models), selected = names(interval_models)[1])
    
  }, ignoreInit = TRUE)
  
  
  observeEvent({
    input$interval_inter_mode
    input$interval_index
    }, {

    # browser()
    
    interaction_mode <- input$interval_inter_mode
    interval_index <- index$custom[[input$interval_index]]

    if(interaction_mode == "Create") {
      
      interval_tmp_mdl$id <- 1
      # update_interval_methods(update_interval_methods() + 1)
      
    } else {}
  }, ignoreInit = TRUE)
  
  
  
  intervalPlot_observer <- observeEvent(event_data("plotly_click", source = "intervalPlot"), {
    
    if (length(interval_tmp_mdl$id) != 0 && input$interval_inter_mode == "Create") {

      # browser()

      click_info <- event_data("plotly_click", source = "intervalPlot")

      if(!is.null(click_info) && click_info[1, "curveNumber"] == 0 && click_info[1, "pointNumber"] != 0) {

        new_id <- click_info[1, "pointNumber"] + 1
        current_id <- interval_tmp_mdl$id

        dup_index <- new_id == current_id

        if(any(dup_index)) {
          interval_tmp_mdl$id <- current_id[-which(dup_index)]

        } else {
          interval_tmp_mdl$id <- sort(c(current_id, new_id))
        }

        shinyjs::runjs("Shiny.setInputValue('plotly_click-intervalPlot', 'null');") # otherwise an event isnt raised when clicking twice on the same point

      } else {
        return()
      }
    } else {
      return()
    }
  }, suspended = TRUE)
  
  
  
  observeEvent(interval_tmp_mdl$id, {
    
    # browser()
    
    new_ids <- interval_tmp_mdl$id
    prev_ids <- interval_tmp_mdl$prev_ids

    req(new_ids)

    del_id <- setdiff(prev_ids, new_ids)

    for (id in del_id){
      removeUI(selector = paste0('#', id), immediate = TRUE)
      runjs(paste0("Shiny.onInputChange('",paste0("sel_",id),"',null)"))
      runjs(paste0("Shiny.onInputChange('",paste0("opt_1_",id),"',null)"))
      runjs(paste0("Shiny.onInputChange('",paste0("opt_2_",id),"',null)"))
    }

    add_id <- sort(setdiff(new_ids, prev_ids))
    
    for (id in add_id) {
      
      if (id == 1) {
        pos <- "placeholder"
      } else {
        pos <- new_ids[which(id == new_ids) - 1]
      }
      
      insertUI(
        selector = paste0("#", pos),
        where = "afterEnd",
        ui = tags$div(
          div(
            div(style="display:inline-block;vertical-align:top; width: 150px;",
                selectInput(paste0("sel_", id), paste0("Model", id),
                            choices = c("None"),
                            selected = "None")),
            div(style="display:inline-block;vertical-align:top; width: 150px;",
                selectInput(paste0("opt_1_", id), "Option 1",
                            choices = c("None"), selected = "None")),
            div(style="display:inline-block;vertical-align:top; width: 150px;",
                selectInput(paste0("opt_2_", id), "Option 2",
                            choices = c("None"), selected = "None"))
          ), id = id), immediate = TRUE
      )
    }

    interval_tmp_mdl$prev_ids <- new_ids

  }, ignoreInit = TRUE)
  
  observe({
    
    req(length(interval_tmp_mdl$id) > 0, !is.null(input$interval_index))
    
    # browser()
    
    method_changed <- FALSE
    
    method_select_choices <- character()
    opt1_select_choices <- list()
    opt2_select_choices <- list()
    
    for (i in seq(length(interval_tmp_mdl$id))) {
      
      id <- interval_tmp_mdl$id[i]
      
      selected_method <- input[[paste0("sel_", id)]]
      selected_opt1 <- input[[paste0("opt_1_", id)]]
      selected_opt2 <- input[[paste0("opt_2_", id)]]
      
      isolate({
        
        sub_index <- interval_tmp_mdl$sub_index()[[i]]
          
        method_select_choices <- "None"
        
        opt1_select_choices[["None"]] <- "None"
        opt2_select_choices[["None"]] <- "None"
      
        if (length(sub_index) == 1) {
          
          method_select_choices <- c(method_select_choices, "Value")
          
          opt1_select_choices[["Value"]] <- "None"
          opt2_select_choices[["Value"]] <- "None"
          
        } else if (length(sub_index) > 1) {
          
          method_select_choices <- c(method_select_choices, "Average",
                                     "Linear Regression", "Blockwise")
          
          point_nb <- length(which(sub_index))
          
          opt1_select_choices[["Average"]] <- "None"
          opt2_select_choices[["Average"]] <- "None"
          
          opt1_select_choices[["Linear Regression"]] <- seq(point_nb - 1)
          opt2_select_choices[["Linear Regression"]] <- c("Force t0", "Force 0")
          
          opt1_select_choices[["Blockwise"]] <- c("Average", "Previous")
          opt2_select_choices[["Blockwise"]] <- "None"
        }
        
        if (is.null(selected_method) || is.null(selected_opt1) || is.null(selected_opt2) || 
            !any(selected_method == method_select_choices)) {
          
          # this is a new or impossible to render select input method
          
          selected_method <- "None"
          
          updateSelectInput(session, paste0("sel_", id),
                            choices = method_select_choices, selected = selected_method)
          updateSelectInput(session, paste0("opt_1_", id),
                            choices = opt1_select_choices[[selected_method]],
                            selected = opt1_select_choices[[selected_method]][1])
          updateSelectInput(session, paste0("opt_2_", id),
                            choices = opt2_select_choices[[selected_method]],
                            selected = opt2_select_choices[[selected_method]][1])
          
          expected_interval_select_input[["method"]][id] <<- "None"
          expected_interval_select_input[["opt1"]][id] <<- "None"
          expected_interval_select_input[["opt2"]][id] <<- "None"
          
          method_changed <- TRUE
          
        } else {
          
          # the method / option select inputs have valid values
          
          if (expected_interval_select_input[["method"]][id] != selected_method ||
              !any(selected_opt1 == opt1_select_choices[[selected_method]]) ||
              !any(selected_opt2 == opt2_select_choices[[selected_method]])
              ) {
            
            # there has been a change in the method select input or selected options are not anymore valid
            
            updateSelectInput(session, paste0("sel_", id),
                              choices = method_select_choices, selected = selected_method)
            updateSelectInput(session, paste0("opt_1_", id),
                              choices = opt1_select_choices[[selected_method]],
                              selected = opt1_select_choices[[selected_method]][1])
            updateSelectInput(session, paste0("opt_2_", id),
                              choices = opt2_select_choices[[selected_method]],
                              selected = opt2_select_choices[[selected_method]][1])
            
            expected_interval_select_input[["method"]][id] <<- selected_method
            expected_interval_select_input[["opt1"]][id] <<- opt1_select_choices[[selected_method]][1]
            expected_interval_select_input[["opt2"]][id] <<- opt2_select_choices[[selected_method]][1]
            
            method_changed <- TRUE
            
          }
          
          else if (expected_interval_select_input[["opt1"]][id] != selected_opt1) {
            
            # there has been a change in the opt1 select input
            
            updateSelectInput(session, paste0("sel_", id),
                              choices = method_select_choices, selected = selected_method)
            updateSelectInput(session, paste0("opt_1_", id),
                              choices = opt1_select_choices[[selected_method]],
                              selected = selected_opt1)
            updateSelectInput(session, paste0("opt_2_", id),
                              choices = opt2_select_choices[[selected_method]],
                              selected = opt2_select_choices[[selected_method]][1])
            
            expected_interval_select_input[["method"]][id] <<- selected_method
            expected_interval_select_input[["opt1"]][id] <<- selected_opt1
            expected_interval_select_input[["opt2"]][id] <<- opt2_select_choices[[selected_method]][1]
            
            method_changed <- TRUE
            
          }
          
          else if (expected_interval_select_input[["opt2"]][id] != selected_opt2) {
            
            # there has been a change in the opt2 select input
            
            updateSelectInput(session, paste0("sel_", id),
                              choices = method_select_choices, selected = selected_method)
            updateSelectInput(session, paste0("opt_1_", id),
                              choices = opt1_select_choices[[selected_method]],
                              selected = selected_opt1)
            updateSelectInput(session, paste0("opt_2_", id),
                              choices = opt2_select_choices[[selected_method]],
                              selected = selected_opt2)
            
            expected_interval_select_input[["method"]][id] <<- selected_method
            expected_interval_select_input[["opt1"]][id] <<- selected_opt1
            expected_interval_select_input[["opt2"]][id] <<- selected_opt2
            
            method_changed <- TRUE
            
          }
          
          else {
            
            # there hasn't been a change in the select inputs, still update in case the choices changed
            
            updateSelectInput(session, paste0("sel_", id),
                              choices = method_select_choices, selected = selected_method)
            updateSelectInput(session, paste0("opt_1_", id),
                              choices = opt1_select_choices[[selected_method]],
                              selected = selected_opt1)
            updateSelectInput(session, paste0("opt_2_", id),
                              choices = opt2_select_choices[[selected_method]],
                              selected = selected_opt2)
            
          }
        }
      })
    }
  })
  
  output$intervalPlot <- renderPlotly({
    
    req(!is.null(interval_tmp_mdl$id))
    req(!is.null(input$interval_index))
    req(!is.null(input$intervalElement))
    
    # browser()
    
    element <- input$intervalElement
    
    plot_time <- convert_time(numeric_time = parsed_info$time_col,
                              to_time_unit = switch(input$timeDisplay,
                                                    Absolute = "Date",
                                                    Relative = input$timeUnit),
                              numeric_origin = parsed_info$time_col[1])
    
    model_time <- parsed_info$delta_t
    
    if (input$interval_inter_mode == "Create") {
      
      index <- index$custom[[input$interval_index]]
      
      plot_x_index <- plot_time[index]
      y_index_value <- analyte()[["value"]][index , element]
      y_index_sd <- analyte()[["SD"]][index , element]
      
      model_x_index <- model_time[index]
      
      plot_t0s <- convert_time(numeric_time = interval_tmp_mdl$slices_t0(),
                               to_time_unit = switch(input$timeDisplay,
                                                     Absolute = "Date",
                                                     Relative = input$timeUnit),
                               numeric_origin = parsed_info$time_col[1])
      model_t0s <- convert_time(numeric_time = interval_tmp_mdl$slices_t0(),
                               to_time_unit = "Seconds",
                               numeric_origin = parsed_info$time_col[1])
      
      plot_tfs <- convert_time(numeric_time = interval_tmp_mdl$slices_tf(),
                               to_time_unit = switch(input$timeDisplay,
                                                     Absolute = "Date",
                                                     Relative = input$timeUnit),
                               numeric_origin = parsed_info$time_col[1])
      model_tfs <- convert_time(numeric_time = interval_tmp_mdl$slices_tf(),
                                to_time_unit = "Seconds",
                                numeric_origin = parsed_info$time_col[1])
        
    }
    else if (input$interval_inter_mode == "View") {
      if (!is.null(interval_models[[input$interval_names]])) {
      }
    }
      intervalPlot <- plot_ly(x= ~plot_x_index, y= ~y_index_value, type = 'scatter', mode = 'markers',
                              source = "intervalPlot", marker = list(size = 10))
      
        if(length(plot_t0s) != 0) {
  
          for (i in seq(length(plot_t0s))) {
              
            id <- interval_tmp_mdl$id[i]
            
            sel_method <- input[[paste0("sel_", id)]]
            sel_opt1 <- input[[paste0("opt_1_", id)]]
            sel_opt2 <- input[[paste0("opt_2_", id)]]
            
            req(!is.null(sel_method), !is.null(sel_opt1), !is.null(sel_opt2))
            req(expected_interval_select_input[["method"]][id] == sel_method,
                expected_interval_select_input[["opt1"]][id] == sel_opt1,
                expected_interval_select_input[["opt2"]][id] == sel_opt2)
            
              plot_t0 <- plot_t0s[i]
              plot_tf <- plot_tfs[i]
              model_t0 <- model_t0s[i]
              model_tf <- model_tfs[i]
              
              sub_index <- interval_tmp_mdl$sub_index()[[i]]
              slice_index <- (model_t0 <= model_time & model_tf >= model_time)
              model_x_slice <- model_time[slice_index]
              plot_x_slice <- plot_time[slice_index]
              
              modeled_data <- get_modeled_data(model_x_index[sub_index],
                                               list(value = y_index_value[sub_index], SD = y_index_sd[sub_index]),
                                               sel_method,
                                               sel_opt1,
                                               sel_opt2,
                                               model_x_slice)
              
              intervalPlot <- intervalPlot %>%
                add_trace(x = rep(plot_t0, times = 2), mode = 'lines+markers',
                          y = c(min(analyte()[["value"]][, element]), max(analyte()[["value"]][, element])),
                          marker = list(size = 1, color = 'rgba(255, 255, 255, 0)')) %>%
                add_trace(x = plot_x_slice,
                          y = modeled_data[["value"]], mode = 'markers',
                          marker = list(size = 8))
          }
        }
      
      intervalPlot_observer$resume()

      intervalPlot

  })
  
  outputOptions(output, "intervalPlot", priority = -1)
  
  # ISTD settings -------------------------------------------------------------
  
  # Display ISTD table
  output$ISTDtable <- DT::renderDT(datatable(
    {

    if(!is.matrix(data_tables$ISTD())) return()

    format(data_tables$ISTD(), digits = 3, scientific=T)

    },
    extensions = c('Scroller', 'FixedColumns', 'Buttons'),
    options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                   scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE, fixedColumns = TRUE,
                   buttons = c('copy', 'csv'),
                   columnDefs = list(list(width = '120px', targets = "_all")))))
  
  # Defines a proxy for changing selections in the table
  ISTDtableProxy <- DT::dataTableProxy("ISTDtable", session = session)
  
  # Keep selection after table update
  observeEvent({
    input$indexISTDchoiceWhat
    input$indexISTDchoiceIn
    input$ISTDinteractionMode
    input$ISTDcorrectionMethod
    }, {
    if (!is.null(input$indexISTDchoiceWhat) && is_valid_ISTD()){
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
  
  # Display blank table
  output$blkTable <- DT::renderDT(datatable({
    
    if(!is.matrix(data_tables$blank())) return()
    
    format(data_tables$blank(), digits = 3, scientific=T)
    
  },
  extensions = c('Scroller', 'FixedColumns', 'Buttons'),
  options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                 scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE, fixedColumns = TRUE,
                 buttons = c('copy', 'csv'),
                 columnDefs = list(list(width = '120px', targets = "_all")))))
  
  # Display blank plot
  output$blankPlot <- renderPlotly({
    
    element <- input$blkPlotElement
    index <- index$custom[[input$blkPlotIndex]]
    
    if(!(element %in% colnames(data_tables$blank()))) return()
    
    x <- make.unique(parsed_info$name_col[index], sep = ".")
    y <- as.numeric(data_tables$blank()[index, element])
    
    blankPlot <- plot_ly(x=x, y=y, type = 'scatter', mode = 'markers')
    
    layout(blankPlot, xaxis = list(type = "category", categoryorder="array", categoryarray=x))
    
  })
  
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
    calibration_parameters <- calibration_model[["modelParam"]]

    calibration_predict <- calibration_data[ ,"concentration"] * calibration_parameters["slope"] +
                           calibration_parameters["intercept"]
    
    signal_residuals <- calibration_data[ ,"value"] - calibration_predict
    
    calibrationPlot <- plot_ly()
    
    if (is.null(calibration_parameters)) {
      
      calibrationPlot <- add_trace(calibrationPlot,
                                   x=NULL, y=NULL, type = 'scatter', mode = 'markers')
      
      calibrationPlot <- add_annotations(calibrationPlot, x = 0.1, y = 0.5,
                                   xref = "paper", yref = "paper",
                                   text = paste("Linear model impossible to compute. Verify that 1/SD or
                                                1/Y weights are not infinite (Y or SD = 0)."),
                                   xanchor = 'left', showarrow = FALSE, font = list(color = "red"))
      
    } else {
      
      calibrationPlot <- add_trace(calibrationPlot,
                                   x=calibration_data[ ,"concentration"],
                                   y=calibration_predict,
                                   type = 'scatter', mode = 'lines')
      
      calibrationPlot <- add_trace(calibrationPlot,
                                   x=calibration_data[ ,"concentration"],
                                   y=calibration_data[ ,"value"],
                                   type = 'scatter', mode = 'markers',
                                   error_y = list(array=calibration_data[ ,"SD"], color = '#000000'))
      
      residualPlot <- plot_ly()
      
      residualPlot <- add_trace(residualPlot,
                                x=calibration_data[,"concentration"],
                                y=signal_residuals,
                                type = 'scatter', mode = 'markers')
      
      subplot(calibrationPlot, residualPlot, nrows = 2)
      
    }
  })
  
  # Drift settings ------------------------------------------------------------
  
  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      set_drift(set_drift() + 1)
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  # Drift index table
  
  output$smpBlkCorTable <- renderTable({
    
    req(index$drift, is.integer(input$e_ind_drift))
    
    drift_ncol_display <- 5
    
    driftSignal <- data_modified$blankCorrectedRatio()[["value"]]
    
    if (input$e_ind_drift < drift_ncol_display){
      low_col <- 1
    }
    else{
      #low_col <- max(which((seq(input$e_ind_drift)) %% drift_ncol_display  == 0))
      low_col <- as.integer((input$e_ind_drift - 1) / drift_ncol_display) * drift_ncol_display + 1
    }
    up_col <- min(low_col + drift_ncol_display - 1, analyteNumber())
    driftTable <- driftSignal[index$drift, low_col:up_col, drop=FALSE]
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
    
    req(index$drift, is.integer(input$e_ind_drift))
    
    driftTime = parsed_info$delta_t[index$drift]
    driftValue = driftData()[["value"]][ ,input$e_ind_drift]
    driftSD = driftData()[["SD"]][ ,input$e_ind_drift]

    driftPlot <- plot_ly()
    driftPlot <- add_trace(driftPlot, x=driftTime, y=driftValue,
                           type = 'scatter', mode = 'markers',
                           error_y = list(array=driftSD, color = '#000000'))
    
    if (process_parameters$driftCorrectedElement[input$e_ind_drift] != "None") {
      
      iTime <- driftTime[1]
      fTime <- tail(driftTime, 1)
      
      timeInterval <- seq(from = as.numeric(iTime), to = as.numeric(fTime),
                          by=as.numeric((fTime - iTime)/100))
      
      driftModel <- data_modified$driftModels()[[input$e_ind_drift]]
      
      driftPredict <- predict(driftModel, newdata = data.frame(driftTime = timeInterval), interval="predict")
      
      driftPlot <- add_trace(driftPlot, x = timeInterval, y = driftPredict[ ,"fit"],
                             type = 'scatter', mode = 'lines')
      
      driftPlot <- add_annotations(driftPlot, x = 0, y = 1,
                                   xref = "paper", yref = "paper",
                                   text = paste("Adjusted R squared: ", summary(driftModel)$adj.r.squared), showarrow = F)
    }
    
    driftPlot
    
  })
  
  observeEvent(input$setDriftCorrection, {
    
    req(index$drift, process_parameters$driftCorrectedElement)

    process_parameters$driftCorrectedElement[input$e_ind_drift] <- input$driftModelSelection
    
  })

  # Data processing -----------------------------------------------------------
  
  # Concentration table
  
  output$conc <- renderTable({
    
    req(data_modified$concentration())
    
    data_modified$concentration()
  })
  
  # Download table
  
  output$downloadData <- downloadHandler(
    filename = paste("data_", input$viewConcentrationIndex, ".csv", sep = ""),
    content = function(file) {
      
      selectedIndex <- index$custom[[input$viewConcentrationIndex]]
      selectedMode <- strtoi(input$viewConcentrationSwitch)
      
      combinedConcSD <- mergeMatrixes(matrix1 = data_modified$concentration()[["value"]],
                                      matrix2 = data_modified$concentration()[["SD"]],
                                      name1=NULL,
                                      name2="SD")

      if (selectedMode <= 2){
        write.csv(data_modified$concentration()[[selectedMode]][selectedIndex, ],
                  file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
      else if (selectedMode == 3){
        write.csv(combinedConcSD[selectedIndex,],
                  file, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    })
  
  
  # Observers -----------------------------------------------------------------
  
  observe({
    if (!is.null(names(index$custom))) {
      
      updateMultipleSelectInput(session = session, input_names = select_input_names,
                                choices = names(index$custom), selected = names(index$custom)[1])
    }
  })
  
  observe({
    updateSelectInput(session,"blkPlotElement",choices=analyteNames(),selected=analyteNames()[1])
    updateSelectInput(session,"calibrationElement",choices=analyteNames(),selected=analyteNames()[1])
    updateSelectInput(session,"e_drift",choices=analyteNames(),selected=analyteNames()[1])
    updateSelectInput(session,"intervalElement",choices=analyteNames(),selected=analyteNames()[1])
  })
  
  # testObs <- observe({
  #   
  #   browser()
  #   
  #   req(extracted() != 0)
  # 
  #   for (id in interval_tmp_mdl$id) {
  # 
  #     i <- which(id == interval_tmp_mdl$id)
  # 
  #     method <- input[[paste0("sel_", id)]]
  #     
  #     isolate({
  #       if (!is.null(method)) {
  #         if (method == "None" || method == "Value" || method == "Average") {
  #           updateSelectInput(session, paste0("opt_1_", id),
  #                             choices = "None", selected = "None")
  #           updateSelectInput(session, paste0("opt_2_", id),
  #                             choices = "None", selected = "None")
  #         } else if (method == "Blockwise") {
  #           
  #           choices <- c("Average", "Previous")
  #           
  #           if (input[[paste0("opt_1_", id)]] %in% choices) {
  #             selected <- input[[paste0("opt_1_", id)]]
  #           } else {
  #             selected <- choices[1]
  #           }
  #           
  #           updateSelectInput(session, paste0("opt_1_", id),
  #                             choices = choices, selected = selected)
  #           updateSelectInput(session, paste0("opt_2_", id),
  #                             choices = "None", selected = "None")
  #         } else if (method == "Linear Regression") {
  #           
  #           point_nb <- length(which(interval_tmp_mdl$sub_index()[[i]]))
  #           
  #           choices_1 <- seq(point_nb - 1)
  #           choices_2 <- c("Force t0", "Force 0")
  #           
  #           if (input[[paste0("opt_1_", id)]] %in% choices_1) {
  #             selected_1 <- input[[paste0("opt_1_", id)]]
  #           } else {
  #             selected_1 <- choices_1[1]
  #           }
  #           
  #           if (input[[paste0("opt_2_", id)]] %in% choices_2) {
  #             selected_2 <- input[[paste0("opt_2_", id)]]
  #           } else {
  #             selected_2 <- choices_2[1]
  #           }
  #           
  #           updateSelectInput(session, paste0("opt_1_", id),
  #                             choices = choices_1, selected = selected_1)
  #           updateSelectInput(session, paste0("opt_2_", id),
  #                             choices = choices_2, selected = selected_2)
  #         }
  #       }
  #     })
  #   }
  #   
  #   testObs$suspend()
  # 
  # }, suspended = TRUE)
  
  # observeEvent(interval_tmp_mdl$interval_sel_inputs()[["method"]], {
  #   
  #   browser()
  #   
  #   for (id in interval_tmp_mdl$id) {
  #     
  #     i <- which(id == interval_tmp_mdl$id)
  #     
  #     method <- input[[paste0("sel_", id)]]
  #     
  #     if (!is.null(method)) {
  #       if (method == "None" || method == "Value" || method == "Average") {
  #         updateSelectInput(session, paste0("opt_1_", id),
  #                           choices = "None", selected = "None")
  #         updateSelectInput(session, paste0("opt_2_", id),
  #                           choices = "None", selected = "None")
  #       } else if (method == "Blockwise") {
  #         
  #         choices <- c("Average", "Previous")
  #         
  #         if (input[[paste0("opt_1_", id)]] %in% choices) {
  #           selected <- input[[paste0("opt_1_", id)]]
  #         } else {
  #           selected <- choices[1]
  #         }
  #         
  #         updateSelectInput(session, paste0("opt_1_", id),
  #                           choices = choices, selected = selected)
  #         updateSelectInput(session, paste0("opt_2_", id),
  #                           choices = "None", selected = "None")
  #       } else if (method == "Linear Regression") {
  #         
  #         point_nb <- length(which(interval_tmp_mdl$sub_index()[[i]]))
  #           
  #         choices_1 <- seq(point_nb - 1)
  #         choices_2 <- c("Force t0", "Force 0")
  #         
  #         if (input[[paste0("opt_1_", id)]] %in% choices_1) {
  #           selected_1 <- input[[paste0("opt_1_", id)]]
  #         } else {
  #           selected_1 <- choices_1[1]
  #         }
  #         
  #         if (input[[paste0("opt_2_", id)]] %in% choices_2) {
  #           selected_2 <- input[[paste0("opt_2_", id)]]
  #         } else {
  #           selected_2 <- choices_2[1]
  #         }
  #         
  #         updateSelectInput(session, paste0("opt_1_", id),
  #                           choices = choices_1, selected = selected_1)
  #         updateSelectInput(session, paste0("opt_2_", id),
  #                           choices = choices_2, selected = selected_2)
  #       }
  #     }
  #   }
  #   
  # }, ignoreInit = TRUE)
  
  observeEvent(input$useBlankCorrection,{
    if(input$useBlankCorrection == TRUE) {
      shinyjs::enable("blkInteractionMode")
    } else {
      shinyjs::disable("blkInteractionMode")
    }
  }, ignoreInit = FALSE)
  
  observe({
    if(input$blkInteractionMode == "process") {
      shinyjs::enable("setBlkInterpolationMethod")
      shinyjs::enable("blkColSlider")
      shinyjs::enable("blkInterpolationMethod")
      shinyjs::enable("indexBlkchoiceIn")
      shinyjs::enable("indexBlkchoiceWhat")
    }
    else if(input$blkInteractionMode == "view") {
      shinyjs::disable("setBlkInterpolationMethod")
      shinyjs::disable("blkColSlider")
      shinyjs::disable("blkInterpolationMethod")
      shinyjs::disable("indexBlkchoiceIn")
      shinyjs::disable("indexBlkchoiceWhat")
    }
  })
  
  observe({
    if(input$ISTDinteractionMode == "process") {
      shinyjs::enable("setISTDcorrectionMethod")
      shinyjs::enable("ISTDcolSlider")
      shinyjs::enable("ISTDcorrectionMethod")
      shinyjs::enable("indexISTDchoiceIn")
      shinyjs::enable("indexISTDchoiceWhat")
    }
    else if(input$ISTDinteractionMode == "view") {
      shinyjs::disable("setISTDcorrectionMethod")
      shinyjs::disable("ISTDcolSlider")
      shinyjs::disable("ISTDcorrectionMethod")
      shinyjs::disable("indexISTDchoiceIn")
      shinyjs::disable("indexISTDchoiceWhat")
    }
  })
}
