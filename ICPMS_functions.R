library(shiny)
library(stringr)
library(matlib)

update_interval_select_inputs <- function() {

  browser()

  for (i in seq(length(interval_tmp_mdl$id))) {

    id <- interval_tmp_mdl$id[i]

    method_select_choices <- "None"

    if (length(interval_tmp_mdl$sub_index()[[i]]) == 1) {
      method_select_choices <- c(method_select_choices, "Value")
    } else if (length(interval_tmp_mdl$sub_index()[[i]]) > 1) {
      method_select_choices <- c(method_select_choices, "Average",
                                 "Linear Regression", "Blockwise")
    }

    selected_method <- input[[paste0("sel_", id)]]

    if (is.null(selected_method) || !any(selected_method == method_select_choices)) {selected_method <- "None"}

    updateSelectInput(session, paste0("sel_", interval_tmp_mdl$id[i]),
                      choices = method_select_choices, selected = selected_method)
    expected_interval_select_input[["method"]][i] <<- selected_method

    if (selected_method == "None" || selected_method == "Value" || selected_method == "Average") {
      updateSelectInput(session, paste0("opt_1_", id),
                        choices = "None", selected = "None")
      updateSelectInput(session, paste0("opt_2_", id),
                        choices = "None", selected = "None")

      expected_interval_select_input[["opt1"]][i] <<- "None"
      expected_interval_select_input[["opt2"]][i] <<- "None"

    } else if (selected_method == "Blockwise") {

      choices <- c("Average", "Previous")

      if (input[[paste0("opt_1_", id)]] %in% choices) {
        selected <- input[[paste0("opt_1_", id)]]
      } else {
        selected <- choices[1]
      }

      updateSelectInput(session, paste0("opt_1_", id),
                        choices = choices, selected = selected)
      updateSelectInput(session, paste0("opt_2_", id),
                        choices = "None", selected = "None")

      expected_interval_select_input[["opt1"]][i] <<- selected
      expected_interval_select_input[["opt2"]][i] <<- "None"

    } else if (selected_method == "Linear Regression") {

      point_nb <- length(which(interval_tmp_mdl$sub_index()[[i]]))

      choices_1 <- seq(point_nb - 1)
      choices_2 <- c("Force t0", "Force 0")

      if (input[[paste0("opt_1_", id)]] %in% choices_1) {
        selected_1 <- input[[paste0("opt_1_", id)]]
      } else {
        selected_1 <- choices_1[1]
      }

      if (input[[paste0("opt_2_", id)]] %in% choices_2) {
        selected_2 <- input[[paste0("opt_2_", id)]]
      } else {
        selected_2 <- choices_2[1]
      }

      updateSelectInput(session, paste0("opt_1_", id),
                        choices = choices_1, selected = selected_1)
      updateSelectInput(session, paste0("opt_2_", id),
                        choices = choices_2, selected = selected_2)

      expected_interval_select_input[["opt1"]][i] <<- selected_1
      expected_interval_select_input[["opt2"]][i] <<- selected_2

    }
  }
}

convert_time <- function(numeric_time, to_time_unit = c("Date", "Seconds", "Minutes", "Hours"),
                         numeric_origin = NULL) {
  
  # numeric_time and numeric_origin must be numeric values equal to the time in seconds elapsed since 1970-01-01 GMT
  
  if (!is.numeric(numeric_time)) {
    stop("numeric_time must be numeric")
  }
  
  if (to_time_unit != "Date" && (is.null(numeric_origin) || !is.numeric(numeric_origin))) {
    stop("numeric_origin must be numeric when to_time_unit != 'Date'")
  }
  
  if (to_time_unit == "Date") {
    date_time <- as.POSIXct(numeric_time, origin = "1970-01-01")
    return(date_time)
  }
  else if (to_time_unit == "Seconds") {
    diff_time <- numeric_time - numeric_origin
  } 
  else if (to_time_unit == "Minutes") {
    diff_time <- (numeric_time - numeric_origin) / 60
  } 
  else if (to_time_unit == "Hours") {
    diff_time <- (numeric_time - numeric_origin) / 3600
  }
  
  return(diff_time)
}

get_modeled_data <- function(x, y, method, opt1, opt2, x_pred) {
  
  # browser()
  
  y_value <- y[["value"]]
  y_sd <- y[["SD"]]
  
  sample_nb <- length(x_pred)
  
  y_pred_value <- rep(NA, times = sample_nb)
  y_pred_sd <- rep(NA, times = sample_nb)
  
  if (method == "None" || is.null(method) || is.null(opt1) || is.null(opt2) ) {
    # do nothing
  }
  else if (method == "Value") {
    
    y_pred_value[1:sample_nb] <- y_value
    y_pred_sd[1:sample_nb] <- y_sd
  }
  else if (method == "Blockwise") {
    
    req(opt1 %in% c("Previous", "Average"), opt2 == "None")
    
    if (opt1 == "Previous") {
      for (i in seq(length(x_pred))) {
        y_pred_value[i] <- y_value[max(which(x_pred[i] >= x))]
        y_pred_sd[i] <- y_sd[max(which(x_pred[i] >= x))]
      }
    }
    else if (opt1 == "Average") {
      for (i in seq(length(x_pred))) {
        if (all(x_pred[i] >= x)) {
          y_pred_value[i] <- mean(tail(y_value, n = 2))
          y_pred_sd[i] <- sd(tail(y_value, n = 2))
        } 
        else if (!any(x_pred[i] >= x)) {
          y_pred_value[i] <- mean(y_value[1:2])
          y_pred_sd[i] <- sd(y_value[1:2])
        }
        else {
          prev_y <- y_value[max(which(x_pred[i] >= x))]
          next_y <- y_value[max(which(x_pred[i] >= x)) + 1]
          y_pred_value[i] <- mean(c(prev_y, next_y))
          y_pred_sd[i] <- sd(c(prev_y, next_y))
        }
      }
    }
  }
  else if (method == "Average") {
    y_pred_value[1:sample_nb] <- mean(y_value)
    y_pred_sd[1:sample_nb] <- sd(y_value)
  }
  else if (method == "Linear Regression") {
    
    req(!is.na(as.numeric(opt1)), opt2 %in% c("Force t0", "Force 0"))
    
    lin_params <- linreg_estimate(x, y_value, as.numeric(opt1))
    lin_prediction <- linreg_predict(x_pred, lin_params[["betas"]],
                                     lin_params[["cov_matrix"]], ratio = FALSE)
    y_pred_value <- lin_prediction[["value"]]
    y_pred_sd <- lin_prediction[["SD"]]
  }
  
  return(list(value = y_pred_value, SD = y_pred_sd))
}

linreg_estimate <- function(x, y, order, w = NULL) {
  
  # browser()
  
  req(order < length(x))
  
  parameters_list <- list(betas = NULL, cov_matrix = NULL)
  model <- lm(y ~ poly(x, order, raw=TRUE))
    
  parameters_list[["betas"]] <- coefficients(model)
  parameters_list[["cov_matrix"]] <- vcov(model)
  
  return(parameters_list)
}

linreg_predict <- function(x, betas, cov_matrix, ratio = FALSE) {
  
  # browser()
  
  prediction_value <- rep(NA, length(x))
  prediction_sd <- rep(NA, length(x))
  
  if(ratio == TRUE && betas[1] == 0) stop("Ratio impossible, beta_0 = 0")
  
  data_nb <- length(x)
  beta_nb <- length(betas)
  order <- beta_nb - 1
  
  X <- matrix(nrow = data_nb, ncol = beta_nb)

  linear_regression_str_expr <- "("
  params <- character()
  
  for (i in 0:order) {
    
    current_param <- paste0("beta_", i)
    params <- c(params, current_param)
    assign(current_param, betas[i + 1])
    
    linear_regression_str_expr <- paste0(linear_regression_str_expr, " + ", current_param, " * x^", i)
    
    X[ , i + 1] <- x^i
  }
  
  linear_regression_str_expr <- paste0(linear_regression_str_expr, ")" )
  
  if(ratio == TRUE) linear_regression_str_expr <- paste0(linear_regression_str_expr, " / beta_0" )
  
  linear_regression_expr <- parse(text = linear_regression_str_expr)
  jacobian_expr <- deriv(expr = linear_regression_expr, params)
  jacobian_values <- attr(eval(jacobian_expr), "gradient")

  prediction_value <- eval(linear_regression_expr)
  
  for (i in seq(data_nb)) {
    jacobian_value <- matrix(jacobian_values[i, ], nrow = 1)
    prediction_sd[i] <- sqrt(as.numeric(jacobian_value%*%cov_matrix%*%t(jacobian_value)))
  }
  
  return(list(value = prediction_value, SD = prediction_sd))
}

updateMultipleSelectInput <- function(session, input_names, choices, selected) {
  for(input_name in input_names) {
    updateSelectInput(session, inputId = input_name, choices = choices, selected = selected)
  }
}

get_parser <- function(parser_name) {
  if (parser_name == "agilent") parse_function_agilent
}

parse_function_agilent <- function(main_file, std_file, ISTD_file) {
  
    header_1 <- scan(main_file, nlines = 1, what = character(),sep=';')
    header_2 <- scan(main_file, nlines = 1, what = character(),sep=';', skip=1)
    
    header_1 <- fill_empty_string(header_1)
    
    main_dataframe <- read.csv2(main_file, skip = 2, header = FALSE, stringsAsFactors = FALSE)
    names(main_dataframe) <- header_1
    
    std_dataframe <- read.table(std_file, header = TRUE, sep=';', stringsAsFactors=FALSE)
    std_elements <- std_dataframe[ , 1]
    
    # Remove duplicated elements
    std_dataframe <- std_dataframe[!duplicated(std_elements), , drop = FALSE]
    
    row.names(std_dataframe) <- std_elements[!duplicated(std_elements)]
    std_dataframe <- std_dataframe[ , -1, drop = FALSE]
    std_matrix <- as.matrix(std_dataframe)
    
    num_col_index <- which(header_2 == "CPS" | header_2 == "CPS RSD")
    
    main_dataframe[ , num_col_index] <- sapply(main_dataframe[ , num_col_index], as.numeric)
    element_value <- as.matrix(main_dataframe[ , header_2 == "CPS", drop = FALSE])
    element_rsd <- as.matrix(main_dataframe[ , header_2 == "CPS RSD", drop = FALSE])
    element_sd <- element_rsd / 100 * element_value
    element_names <- colnames(element_value)
    
    analyte_index <- !grepl("ISTD", colnames(element_value))
    ISTD_index <- !analyte_index
    
    smp_nb <- nrow(main_dataframe)
    time_col <- as.numeric(as.POSIXct(main_dataframe[ , which(header_2 == "Acq. Date-Time")], format="%d/%m/%Y %H:%M"))
    name_col <- main_dataframe[ , which(header_2 == "Sample Name")]
    type_col <- main_dataframe[ , which(header_2 == "Type")]
    lvl_col <- main_dataframe[ , which(header_2 == "Level")]
    delta_t <- time_col - time_col[1]
    
    row.names(element_value) <- name_col
    row.names(element_rsd) <- name_col
    row.names(element_sd) <- name_col
    
    return(list(element_value = element_value, element_sd = element_sd, element_names = element_names,
                std = std_matrix,
                analyte_index = analyte_index, ISTD_index = ISTD_index, num_col_index = num_col_index,
                smp_nb = smp_nb, time_col = time_col, name_col = name_col, type_col = type_col, lvl_col = lvl_col, delta_t = delta_t))
}

agilentElementNamePattern <- "[ ]{2}[A-Z]{1}[a-z]*[ ]{2}"

getConcentration <- function(m, Cm, sig, vsig, dft, vdft) {
  #Function parameters:
  #m and Cm are the calibration linear regression parameters and their covariance matrix
  #sig and vsig are the signal and signal variance (signal can be counts, ratio, can be blank corrected or not)
  #dft and vdft are the drift factor and drift factor variance (the signal is drift corrected using the signal/drift factor ratio)
  #Returns:
  #concentrations, a matrix containing the concentrations of the the samples in the first column, and the standard error in the second column
  
  #m is of class matrix and of dimension (2,1) and contains the calibration linear regression parameters as double
  #m is of class matrix and of dimension (2,2) and contains the variance and covariance of calibration linear regression parameters as double
  
  intercept <- as.numeric(m[1])
  slope <- as.numeric(m[2])
  
  #sig, vsig, dft, vdft are of similar class and size, they are numeric of length sampleNumber
  sampleNumber <- length(sig)
  
  concentrations <- matrix(rep(NA, sampleNumber * 2), nrow = sampleNumber, ncol = 2)
  row.names(concentrations) <- row.names(sig)
  colnames(concentrations) <- c("value", "SD")
  
  #The concentration of each samples is calculated independently with this loop. Results are concatenated in the concentration matrix
  for (i in seq(sampleNumber)) {
    #J is of class matrix and of dimension (1,4) 
    #J is the Jacobian matrix associated with parameters intercept, slope, signal, drift in the equation concentration = (signal/drift - intercept)/slope
    #J is used for the calculation of the standard errors on the concentrations
    J <- t(c(-1/slope,-((sig[i]-intercept)/(dft[i]*(slope^2))),1/(dft[i]*slope),-sig[i]/(dft[i]^2*slope)))
    #Cy is of class matrix and of dimension (4,4) 
    #Cy is the covariance matrix for the parameters intercept, slope, signal, drift
    #Cy is used for the calculation of the standard errors on the concentrations
    Cy <- cbind(Cm,matrix(0,2,2))
    Cy <- rbind(Cy,matrix(0,2,4))
    Cy[3,3] <- vsig[i]
    Cy[4,4] <- vdft[i]
    
    #concentrations is of class matrix and of dimension (sampleNumber,2) at the end of the for loop (addition of one row per loop) 
    #concentrations contains the sample concentrations and their standard error
    concentrations[i , ] <- c((sig[i]/dft[i]-intercept) / slope, sqrt(as.numeric(J%*%Cy%*%t(J))))
  }
  return(concentrations)
}

getCalibrationModel <- function(element, processParameters, calibrationData) {
  #Function parameters:
  #element is the element full name such as it appears in the header of the original csv file
  #processParameters is the reactiveValue that contains the parameters fi, wr and w (force intercept, weighted regression and weights)
  #processParameters$fi and $wr are logical vectors of length analyte number
  #calibrationData contains the result of the getCalibrationData function, a matrix of nrow StdNumber of ncol 3 ("Value", "SD", "Concentration")
  #Returns:
  #A list containing the linear regression model parameters as "modelParam" and the covariance matrix as "covMatrix"
  
  
  
  fi <- processParameters$forceIntercept
  wr <- processParameters$useWeithedRegression
  w <- processParameters$regressionWeight
  
  stdValue <- calibrationData[, "value"]
  stdConcentrations <- calibrationData[, "concentration"]
  if (fi[[element]] == TRUE) {
    dfRegression <- length(stdConcentrations) - 1
    
    intercept <- calibrationData[1,"value"]
    stdValue <- stdValue - intercept
    
    G <- stdConcentrations
  }
  else {
    dfRegression <- length(stdConcentrations) - 2
    
    G <- matrix(cbind(rep(1, length(stdConcentrations)), stdConcentrations), ncol=2)# 1-Intercept 2-Slope
  }
  
  if (wr[[element]] == TRUE) {
    calibrationWeights <- getWeights(calibrationData = calibrationData, fn = w[[element]])
    
    if (Inf %in% calibrationWeights) return(list(modelParam = NULL, covMatrix = NULL))
    
    Cd <- diag(x=calibrationWeights, nrow=length(stdConcentrations), ncol=length(stdConcentrations))
  }
  else {
    Cd <- diag(x=1, nrow=length(stdConcentrations), ncol=length(stdConcentrations))
  }

  m <- solve(t(G)%*%solve(Cd)%*%G)%*%t(G)%*%solve(Cd)%*%stdValue
  chi2 <- as.vector((t(stdValue-G%*%m)%*%solve(Cd)%*%(stdValue-G%*%m)))/dfRegression
  Cm <- solve(t(G)%*%solve(Cd)%*%G)*chi2
  
  if (fi[[element]] == TRUE) {
    m <- c(intercept,m)
    Cslope <- Cm
    Cm <- matrix(0,2,2)
    Cm[2,2] <- Cslope
  }
  
  names(m) <- c("intercept","slope")
  
  return(list(modelParam = m, covMatrix = Cm))
}

subsetDfList <- function(dfList, rowLogical, columnLogical) {
  subDfList <- list()
  for (i in 1:length(dfList)) {
    subDfList[[i]] <- dfList[[i]][rowLogical, columnLogical]
    names(subDfList)[i] <- names(dfList)[i]
  }
  return(subDfList)
}

getUncertaintyInterval <- function(signal) {
  return(list(uBound=signal[["value"]] + signal[["SD"]], lBound=signal[["value"]] - signal[["SD"]]))
}

getModifierName <- function(modifier) {
    return(paste(modifier[["method"]], " lines:", modifier[["whatIndex"]][1], ":", tail(modifier[["whatIndex"]], n=1), sep=""))
}

getModifiedData <- function(dat, modifiers, selectedModifiers) {
  if (length(modifiers) == 0) {return(dat)}
  else {
    
    dat_mod <- dat
    for (s in selectedModifiers) {
      mod <- modifiers[[s]]
      newValues <- replaceValues(dat[["value"]], dat[["SD"]], mod[["elements"]], mod[["method"]], mod[["whatIndex"]], mod[["inIndex"]])
      dat_mod[["value"]][mod[["whatIndex"]],] <- newValues[["value"]][mod[["whatIndex"]],]
      dat_mod[["SD"]][mod[["whatIndex"]],] <- newValues[["SD"]][mod[["whatIndex"]],]
    }
    return(dat_mod)
  }
}

getWeights <- function(calibrationData, fn) {
  if (fn == "1/var") {
    return(1/(calibrationData[,"SD"])^2)
  }
  else if (fn == "1/Y") {
    return(1/(calibrationData[,"value"]))
  }
  else if (fn == "1/max(var,Y)") {
    bindedCols <- cbind((calibrationData[,"SD"])^2, calibrationData[,"value"])
    maxCol <- apply(X=bindedCols,MARGIN=1,FUN=max)
    return(1/(maxCol))
  }
  else {
    return()
  }
}

getElementDriftIndex <- function(elementFullName, stdDataFrame, stdIdentificationColumn, driftIndex) {
  
  #Extracts the element name from the element full name
  eName <- getElementName(elementFullName, agilentElementNamePattern)
  
  #Rises an error if the element name is not found in the standard dataframe
  if (!(eName %in% row.names(stdDataFrame))) {
    stop(paste(eName, " not found in the standard file. If not in the standards, create an empty line for ", eName, " in the file.", sep=""))
  } 
  else {}
  
  stdIndex <- which(!is.na(stdDataFrame[eName,]))
  
  if (length(stdIndex) != 0){
    firstStandardNumIndex <- which(colnames(stdDataFrame)[stdIndex[1]] == make.names(stdIdentificationColumn))
    eDritftIndex <- driftIndex[driftIndex >= firstStandardNumIndex]
    return(eDritftIndex)
  }
  
  else{
    return(NA)
  }
}

getCalibrationData <- function(elementFullName, signal, stdIdentificationColumn, stdDataFrame, truncate = FALSE, smpIndex = NULL) {
  
  #Extracts the element name from the element full name
  eName <- getElementName(elementFullName, agilentElementNamePattern)
  
  #Rises an error if the element name is not found in the standard dataframe
  if (!(eName %in% row.names(stdDataFrame))) {
    stop(paste(eName, " not found in the standard file. If not in the standards, create an empty line for ", eName, " in the file.", sep=""))
  } 
  else {}
  
  #stdConcentration contains the concentrations of the element eName in the standards
  stdConcentration <- t(stdDataFrame[eName, , drop = FALSE])
  #if the element eName has concentrations associated to it, calibrate
  if (sum(stdConcentration, na.rm=TRUE) != 0)
  {
    #Discard non numeric values in stdConcentration
    stdConcentration <- stdConcentration[!is.na(stdConcentration), 1, drop = FALSE]
    
    #Searches for the standard names (row names in stdConcentration) that are associated with numeric values in
    #the level column to create an index that will be used to find the signal associated to the standards
    calibrationStdNumIndex <- sapply(paste("^", row.names(stdConcentration), "$", sep=""),
                                     grep, make.names(stdIdentificationColumn))
    
    #Rises an error if one of the standard id is not found (if grep didn't find a standard name, return a 0 length integer)
    if (!all(sapply(calibrationStdNumIndex,length) > 0)) {
      stop(paste("For ", eName, ", no name match between standard file and standard identifiers in the sequence.", sep=""))
    } 
    else {} 
    
    eValue <- signal[[1]][calibrationStdNumIndex, elementFullName]
    eSD <- signal[[2]][calibrationStdNumIndex, elementFullName]
  }
  else
  {
    return(NA)
  }
  
  return(cbind(value=eValue, SD=eSD, concentration=as.numeric(stdConcentration)))
}

getElementName <- function(elementFullName, pattern) {
  return(gsub(" ", "", str_extract(elementFullName, pattern), fixed = TRUE))
}

mergeMatrixes <- function(matrix1, matrix2, name1=NULL, name2=NULL) {
  if (!all(dim(matrix1) == dim(matrix2))){
    stop("Matrixes of different sizes")
  }
  
  columnNumber <- ncol(matrix1)
  
  if (is.null(name1)){}
  else{colnames(matrix1) <- rep(name1,columnNumber)}
  
  if (is.null(name2)){}
  else{colnames(matrix2) <- rep(name2,columnNumber)}
  
  for(i in seq(columnNumber)){
    if (i == 1){
      combinedMatrix <- cbind(matrix1[,1],matrix2[,1])
      colnames(combinedMatrix) <- c(colnames(matrix1)[1],colnames(matrix2)[1])
    }
    else{
      combinedMatrix <- cbind(combinedMatrix, matrix1[,i], matrix2[,i])
      colnames(combinedMatrix)[ncol(combinedMatrix) - 1] <- colnames(matrix1)[i]
      colnames(combinedMatrix)[ncol(combinedMatrix)] <- colnames(matrix2)[i]
    }
  }
  
  return(combinedMatrix)
}


##Returns the state as a text message. State is TRUE or FALSE.
##If state is TRUE, either stateTxt and validStateTxt is returned, or nothing is returned (used for warnings to pop only when the state is invalid)
##If state is FALSE, stateTxt and invalidStateTxt is returned
renderState <- function(state, stateTxt, invalidStateTxt, validStateTxt) {
  
  if (state == TRUE){
    #This line is used for warnings to pop up only when invalidState
    if (is.null(validStateTxt)){return()}
    else {displayedTxt <- paste('<span style=\"color:green\">', stateTxt, validStateTxt)}}
  else{displayedTxt <- paste('<span style=\"color:red\">', stateTxt, invalidStateTxt)}
  
  return(displayedTxt)
}

search_closest_in_index <- function(search_position, search_in, search_where) {
  
  if (search_where == "previous"){
    if (TRUE %in% (search_in <= search_position)){
      return(max(search_in[search_in <= search_position]))
    } else {
      return(NULL)
    }
  }
  else if (search_where == "next"){
    if (TRUE %in% (search_in >= search_position)){
      return(min(search_in[search_in >= search_position]))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
    }
}

propagateUncertainty <- function(a, b, operation){
  
  a_value = a[[1]]
  b_value = b[[1]]
  a_SD = a[[2]]
  b_SD = b[[2]]
  
  if (operation == "addition"){
    value <- a_value + b_value
    SD <- sqrt((a_SD)^2+(b_SD)^2)
  }
  else if (operation == "substraction"){
    value <- a_value - b_value
    SD <- sqrt((a_SD)^2+(b_SD)^2)
  }
  else if (operation == "multiplication"){
    value <- a_value * b_value
    SD <- sqrt((a_SD/a_value)^2+(b_SD/b_value)^2)*value
  }
  else if (operation == "division"){
    value <- a_value / b_value
    SD <- sqrt((a_SD/a_value)^2+(b_SD/b_value)^2)*value
  }
  else {return(NULL)}
  
  return(list(value = value , SD=  SD))
}

##Function that replaces each value in a text vector by its previous value if the current value is empty (=="")
##Note that if the first value of the vector is empty, there will be empty values in the returned vector
fill_empty_string <- function(char_vector){
  for (i in 2:length(char_vector)){
    if (char_vector[i] == ""){
      char_vector[i] <- char_vector[i-1]
    }
  }
  
  return(char_vector)
}

rm_duplicate_lines <- function(df){
  n_line <- 1
  while (n_line < nrow(df)){
    linesToDelete <- numeric()
    for (i in (n_line + 1):nrow(df)){
      if (identical(as.character(df[i,]), as.character(df[n_line,]))){
        linesToDelete <- c(linesToDelete, i)
      }
      else{}
    }
    
    if (length(linesToDelete) != 0){
      df <- df[-linesToDelete, ]
    }
    else {}
    
    n_line <- n_line + 1
  }
  
  return(df)
}

# Replaces specific values of a matrix based on reference values
replaceValues <- function(value, SD, elements, replace_how, replace_what, replace_in){

  browser()
  
  if (is.null(value) | is.null(SD)) return(NULL)
  
  if (replace_how == "none" | is.null(replace_what) | is.null(replace_in) | is.null(elements)) return(list(value=value,SD=SD))
  
  if (!is.numeric(replace_what) | !is.numeric(replace_in)) {
    stop("In function replace_values, non-numerical index.")
    return(NULL)
  }
  
  line_nb <- nrow(value)
  
  # Replace in selected "All" index case
  if (length(replace_in) == line_nb){
    replace_in <- seq(line_nb)[-replace_what]
  }
  
  # Replaced values cannot be in reference values
  if (any(replace_what %in% replace_in)) return(list(value=value,SD=SD))

  for (element in elements){
    for (n_line in replace_what){
      # Replaces the values in the signal and SD matrix based on replacement type
      previous_line <- search_closest_in_index(search_position = n_line,
                                          search_in = replace_in,
                                          search_where = "previous")
      previous_value <- value[previous_line, element]
      previous_SD <- SD[previous_line, element]

      if (replace_how == "mean"){
        next_line <- search_closest_in_index(n_line, replace_in, "next")
        next_value <- value[next_line, element]
        
        value[n_line, element] <- mean(c(previous_value, next_value))
        SD[n_line, element] <- sd(c(previous_value, next_value))
      } 
      else if (replace_how == "previous"){
        value[n_line, element] <- previous_value
        SD[n_line, element] <- previous_SD
      }
    }
  }
  return(list(value = value, SD = SD))
}

##Function to create the ISTD template based on main datafile
##Returns a dataframe containing a column analyte with the analyte names
##and a column ISTD with the ISTD names
createISTDtemplate <- function(dataFileName){

  #header_1 contains the first line of the file
  header_1 <- scan(dataFileName, nlines = 1, what = character(),sep=';')
  #header_2 contains the second line of the file
  header_2 <- scan(dataFileName, nlines = 1, what = character(),sep=';',skip=1)
  
  #fills empty values in header_1 based on previous values
  header_1 <- fill_empty_string(header_1)
  
  #finds the analyte and the ISTD in the headers based on the occurence of CPS and ISTD keywords
  analyteIndex <- (header_2 == "CPS") & (!grepl("ISTD", header_1))
  ISTDindex <- (header_2 == "CPS") & (grepl("ISTD", header_1))
  
  analyteColumn <- header_1[analyteIndex]
  ISTD <- header_1[ISTDindex]
  ISTDColumn <- c(ISTD, rep(NA, times=length(analyteColumn) - length(ISTD)))

  return(data.frame(analyte=analyteColumn, ISTD=ISTDColumn))
}
  
createISTDMatrix <- function(ISTD_file, ISTDcount){
  
  if(is.null(ISTD_file) | is.null(ISTDcount)){return(1)}
  
  for (j in 1:nrow(ISTD_file)){
    eISTD <- ISTD_file[j, 2]
    if (j == 1){
      ISTDmatrixValue <- ISTDcount[["value"]][ , eISTD, drop = FALSE]
      ISTDmatrixSD <- ISTDcount[["SD"]][ , eISTD, drop = FALSE]
      print(2)
    } else {
      ISTDmatrixValue <- cbind(ISTDmatrixValue, ISTDcount[["value"]][ , eISTD, drop = FALSE])
      ISTDmatrixSD <- cbind(ISTDmatrixSD, ISTDcount[["SD"]][ , eISTD, drop = FALSE])
    }
  }
  return(list(value=ISTDmatrixValue, SD=ISTDmatrixSD))
}