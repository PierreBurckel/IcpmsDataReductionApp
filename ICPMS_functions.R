library(shiny)
library(stringr)
library(matlib)

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
  concentrations <- numeric()
  
  #sig, vsig, dft, vdft are of similar class and size, they are numeric of length sampleNumber
  sampleNumber <- length(sig)
  
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
    concentrations <- rbind(concentrations, t(c(value = (sig[i]/dft[i]-intercept)/slope, SE = sqrt(as.numeric(J%*%Cy%*%t(J))))))
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
    #browser()
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
  stdConcentration <- t(stdDataFrame[eName,])
  #if the element eName has concentrations associated to it, calibrate
  if (sum(stdConcentration, na.rm=TRUE) != 0)
  {
    #Discard non numeric values in stdConcentration
    stdConcentration <- stdConcentration[!is.na(stdConcentration),1, drop = FALSE]
    
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

##Function that converts a numerical index to a logical index
##Takes as argument the numerical index and the number of line expected in the logical index
##Returns the logical index
convertToLogical <- function(nIndex, lineNb){
  lIndex <- rep(FALSE, lineNb)
  lIndex[nIndex] <- TRUE
  
  return(lIndex)
}
##Function that returns the closest number in an array (numericalIndex) before or after a reference position
##If there is no value matching the request, returns NULL
getClosestInIndex <- function(refPosition, numericalIndex, searchWhere){
  
  if (searchWhere == "prev"){
    if (TRUE %in% (numericalIndex <= refPosition)){
      return(max(numericalIndex[numericalIndex <= refPosition]))
    }
    else{
      return(NULL)
    }
  }
  else if (searchWhere == "next"){
    if (TRUE %in% (numericalIndex >= refPosition)){
      return(min(numericalIndex[numericalIndex >= refPosition]))
    }
    else{
      return(NULL)
    }
  }
  else{return(NULL)}
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
  
  return(list(value=as.data.frame(value), SD=as.data.frame(SD)))
}

##Function that replaces each value in a text vector by its previous value if the current value is empty (=="")
##Note that if the first value of the vector is empty, there will be empty values in the returned vector
fillEmptyStrings <- function(textVector){
  for (i in 2:length(textVector)){
    if (textVector[i] == ""){
      textVector[i] <- textVector[i-1]
    }
  }
  
  return(textVector)
}

removeDuplicateLines <- function(df){
  nLine <- 1
  while (nLine < nrow(df)){
    linesToDelete <- numeric()
    for (i in (nLine + 1):nrow(df)){
      if (identical(as.character(df[i,]), as.character(df[nLine,]))){
        linesToDelete <- c(linesToDelete, i)
      }
      else{}
    }
    
    if (length(linesToDelete) != 0){
      df <- df[-linesToDelete, ]
    }
    else {}
    
    nLine <- nLine + 1
  }
  
  return(df)
}
##Function that replaces all values in lineIndex
##Values are replaced by new values, based on the sourceLineIndex index
##If the sourceLineIndex are set to NULL, the lines that are not in the lineIndex (i.e. !lineIndex) are used as source
##The new values can be the mean and SD of previous and following sourceLineIndex or the value and SD of the previous sourceLine
##ICPsignal and SD are the dataframe containing the signal and SD of the signal 
##col_range is the column range for line index replacement,
##replace_type is the replacement type (mean or previous value)
replaceValues <- function(ICPvalue, SD, colRange, replaceType, lineIndex, sourceLineIndex){
  #This condition is useful when the signal hasn't been extracted and the function is called -> returns null
  if(is.null(ICPvalue)| is.null(SD)){return(NULL)}
  #This condition is useful in shiny tables to return the table signal and SD when no replacements are required
  if (replaceType == "none" | is.null(lineIndex) | is.null(colRange)) {return(list(value=ICPvalue,SD=SD))}
  #Stores the line number of the signal whether it is a dataframe (nrow) or a vector (length)
  if(is.integer(nrow(ICPvalue))){
    lineNb <- nrow(ICPvalue)
  } else {
    lineNb <- length(ICPvalue)
  }
  #The function requires a logical index. If the index is non-logical (i.e. numeric) it is converted to a logical index
  if (!is.logical(lineIndex)) {
    lineIndex <- convertToLogical(lineIndex, lineNb)
  }
  #The function accepts a NULL sourceLineIndex. The sourceLineIndex variable takes the opposite value of the logical lineIndex (!lineIndex)
  if (is.null(sourceLineIndex)){
    sourceLineIndex <- !lineIndex
  }
  #As for lineIndex, if sourceLineIndex is numerical it needs to be converted to a logical index
  else if (!is.logical(sourceLineIndex)) {
    sourceLineIndex <- convertToLogical(sourceLineIndex, lineNb)
  }
  #If the source line index is the "all" index, converts the index to the inverse of the lineIndex
  if (all(sourceLineIndex == rep(TRUE, rep=lineNb))){
    sourceLineIndex <- !lineIndex
  }

  #Now that we are sure that lineIndex and sourceLineIndex are logical, we create their numerical equivalent
  nLineIndex <- which(lineIndex)
  nSourceLineIndex <- which(sourceLineIndex)
  #Apply modifications only on colRange columns
  for (j in colRange){
    #Apply modifications only on lines in the lineIndex
    for (nLine in nLineIndex){
      #Replaces the values in the signal and SD matrix based on replacement type
      prevLine <- getClosestInIndex(nLine, nSourceLineIndex, "prev")
      prevValue <- ICPvalue[prevLine,j]
      prevValueSD <- SD[prevLine,j]

      if (replaceType == "mean"){
        nextLine <- getClosestInIndex(nLine, nSourceLineIndex, "next")
        nextValue <- ICPvalue[nextLine,j]
        ICPvalue[nLine,j] <- mean(c(prevValue,nextValue))
        SD[nLine,j] <- sd(c(prevValue,nextValue))
      } 
      else if (replaceType == "prev"){
        ICPvalue[nLine,j] <- prevValue
        SD[nLine,j] <- prevValueSD
      }
    }
  }
  return(list(value=ICPvalue,SD=SD))
}

##Function to create the ISTD template based on raw datafile
##Returns a dataframe containing a column analyte with the analyte names
##and a column ISTD with the ISTD names
createISTDtemplate <- function(dataFileName){

  #header_1 contains the first line of the file
  header_1 <- scan(dataFileName, nlines = 1, what = character(),sep=';')
  #header_2 contains the second line of the file
  header_2 <- scan(dataFileName, nlines = 1, what = character(),sep=';',skip=1)
  
  #fills empty values in header_1 based on previous values
  header_1 <- fillEmptyStrings(header_1)
  
  #finds the analyte and the ISTD in the headers based on the occurence of CPS and ISTD keywords
  analyteIndex <- (header_2 == "CPS") & (!grepl("ISTD", header_1))
  ISTDindex <- (header_2 == "CPS") & (grepl("ISTD", header_1))
  
  analyteColumn <- header_1[analyteIndex]
  ISTD <- header_1[ISTDindex]
  ISTDColumn <- c(ISTD, rep(NA, times=length(analyteColumn) - length(ISTD)))

  return(data.frame(analyte=analyteColumn, ISTD=ISTDColumn))
}

##Function that takes the raw data and standard file along with some options
##and returns the extracted information from these file. Extracted information
##is CPS, SD, IS.CPS, IS.SD, etc...
extractData <- function(dataFileName, stdFileName){
  
  #header_1 contains the first line of the file
  header_1 <- scan(dataFileName, nlines = 1, what = character(),sep=';')
  #header_2 contains the second line of the file
  header_2 <- scan(dataFileName, nlines = 1, what = character(),sep=';', skip=1)
  
  #fills empty values in header_1 based on previous values
  header_1 <- fillEmptyStrings(header_1)
  
  #Import raw data without headers as they are imported previously. header_1 is used for column names
  dat.raw <- read.csv2(dataFileName, skip = 2, header = FALSE, stringsAsFactors=FALSE)
  names(dat.raw) <- header_1
  
  #Import the standard data, we consider here that there are headers corresponding to the standard names
  StdDataframe <- read.table(stdFileName, header = TRUE, sep=';', stringsAsFactors=FALSE)
  
  #Remove potential duplicates of element names in the standard dataframe
  StdDataframe <- removeDuplicateLines(StdDataframe)
  
  #Assigns the first column of the std file to the row names then removes the first column
  row.names(StdDataframe) <- StdDataframe[,1]
  StdDataframe <- StdDataframe[,2:length(StdDataframe)]
  
  return(list(raw=dat.raw, std=StdDataframe, header_1=header_1, header_2=header_2))
}
  
  
createISTDMatrix <- function(ISTD_file, ISTDcount){

  if(is.null(ISTD_file) | is.null(ISTDcount)){return(1)}
  
  for (j in 1:nrow(ISTD_file)){
    eISTD <- ISTD_file[j,2]
    if (j == 1){
      ISTDmatrixValue <- ISTDcount[["value"]][eISTD]
      ISTDmatrixSD <- ISTDcount[["SD"]][eISTD]
      print(2)
    } else {
      ISTDmatrixValue <- cbind(ISTDmatrixValue, ISTDcount[["value"]][eISTD])
      ISTDmatrixSD <- cbind(ISTDmatrixSD, ISTDcount[["SD"]][eISTD])
    }
  }
  return(list(value=ISTDmatrixValue, SD=ISTDmatrixSD))
}