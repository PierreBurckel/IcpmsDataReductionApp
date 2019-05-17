library(shiny)
library(stringr)

agilentElementNamePattern <- "[ ]{2}[A-Z]{1}[a-z]*[ ]{2}"

getWeights <- function(calibrationData, fn) {
  if (fn == "1/SD") {
    return(1/(calibrationData[,"RSD"]/100 *calibrationData[,"Signal"]))
  }
  else if (fn == "1/SD^2") {
    return(1/(calibrationData[,"RSD"]/100 *calibrationData[,"Signal"])^2)
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
    firstStandardNumIndex <- grep(colnames(stdDataFrame)[stdIndex[1]],make.names(stdIdentificationColumn))
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
    
    eSignal <- signal[[1]][calibrationStdNumIndex, elementFullName]
    eRSD <- signal[[2]][calibrationStdNumIndex, elementFullName]
  }
  else
  {
    return(NA)
  }
  
  return(cbind(Signal=eSignal, RSD=eRSD, Concentration=as.numeric(stdConcentration)))
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

##Function that calculates the RSD of a collection of numerical values
calculateRSD <- function(values){
  if (length(values) > 1){
    return(sd(values)/mean(values)*100)
  }
  else{
    return(NULL)
  }
}


propagateUncertainty <- function(a, b, operation){
  
  a_signal = a[[1]]
  b_signal = b[[1]]
  a_RSD = a[[2]]
  b_RSD = b[[2]]
  if (operation == "addition"){
    signal <- a_signal + b_signal
    RSD <- sqrt((a_RSD/100*a_signal)^2+(b_RSD/100*b_signal)^2)/signal*100
  }
  else if (operation == "substraction"){
    signal <- a_signal - b_signal
    RSD <- sqrt((a_RSD/100*a_signal)^2+(b_RSD/100*b_signal)^2)/signal*100
  }
  else if (operation == "multiplication"){
    signal <- a_signal * b_signal
    RSD <- sqrt((a_RSD/100)^2+(b_RSD/100)^2)*100
  }
  else if (operation == "division"){
    signal <- a_signal / b_signal
    RSD <- sqrt((a_RSD/100)^2+(b_RSD/100)^2)*100
  }
  else {return(NULL)}
  
  return(list(signal=signal, RSD=RSD))
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
##ICPsignal and RSD are the dataframe containing the signal and RSD of the signal 
##col_range is the column range for line index replacement,
##replace_type is the replacement type (mean or previous value)
replaceValues <- function(ICPsignal, RSD, colRange, replaceType, lineIndex, sourceLineIndex){
  #This condition is useful when the signal hasn't been extracted and the function is called -> returns null
  if(is.null(ICPsignal)| is.null(RSD)){return(NULL)}
  #This condition is useful in shiny tables to return the table signal and RSD when no replacements are required
  if (replaceType == "none" | is.null(lineIndex) | is.null(colRange)) {return(list(ICPsignal,RSD))}
  #Stores the line number of the signal whether it is a dataframe (nrow) or a vector (length)
  if(is.integer(nrow(ICPsignal))){
    lineNb <- nrow(ICPsignal)
  } else {
    lineNb <- length(ICPsignal)
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
      #Replaces the values in the signal and RSD matrix based on replacement type
      prevLine <- getClosestInIndex(nLine, nSourceLineIndex, "prev")
      prevValue <- ICPsignal[prevLine,j]
      prevValueRSD <- RSD[prevLine,j]

      if (replaceType == "mean"){
        nextLine <- getClosestInIndex(nLine, nSourceLineIndex, "next")
        nextValue <- ICPsignal[nextLine,j]
        ICPsignal[nLine,j] <- mean(c(prevValue,nextValue))
        RSD[nLine,j] <- calculateRSD(c(prevValue,nextValue))
      } 
      else if (replaceType == "prev"){
        ICPsignal[nLine,j] <- prevValue
        RSD[nLine,j] <- prevValueRSD
      }
    }
  }
  return(list(signal=ICPsignal,RSD=RSD))
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
##is CPS, RSD, IS.CPS, IS.RSD, etc...
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
  
  
createISTDMatrix <- function(ISTD_file, ISTD_signal){
  
  if(is.null(ISTD_file) | is.null(ISTD_signal)){return(1)}
  
  for (j in 1:nrow(ISTD_file)){
    eISTD <- ISTD_file[j,2]
    if (j == 1){
      ISTD.matrix <- ISTD_signal[eISTD]
    } else {
      ISTD.matrix <- cbind(ISTD.matrix, ISTD_signal[eISTD])
    }
  }
  return(ISTD.matrix)
}


processData <- function(signal, eFullNames, StdDataframe, drift_ind, levelColumn, timeColumn, eDriftChoice){
  
  elementNumber <- length(eFullNames)
  
  for (i in 1:elementNumber){

    eFullName <- eFullNames[i]

    eCalibrationData <- getCalibrationData(elementFullName = eFullName, signal=signal,
                                       stdIdentificationColumn=levelColumn, stdDataFrame = StdDataframe)

    eDriftIndex <-  getElementDriftIndex(elementFullName = eFullName, stdDataFrame = StdDataframe, 
                                        stdIdentificationColumn=levelColumn, driftIndex = drift_ind)

    dtimeColumn <- timeColumn - timeColumn[eDriftIndex][1]
    dtimeColumn[dtimeColumn < 0] <- 0

    eDriftData <- cbind(Signal=signal[[1]][eDriftIndex,i],dt=dtimeColumn[eDriftIndex])  

    #defines the calibration linear model
    calibrationModel <- lm(Concentration ~ 0+Signal, data=as.data.frame(eCalibrationData))
    
    #Here be drift model creation and signal prediction along the whole sequence for the drift standard
    #if no drift signal, don't try to create a model for the drift and set drift prediction to NA
    #else create a model. Predict the drift signal on the whole sequence. Doesn't mean that it is going to be used
    if (all(is.na(signal[[1]][eDriftIndex,i]))){
      driftPredict = rep(NA, nrow(signal[[1]]))
    }
    else{
      driftModel <- lm(Signal ~ poly(dt, degree=2, raw=TRUE), data = as.data.frame(eDriftData))
      driftPredict=predict(driftModel, newdata = data.frame(dt=dtimeColumn))
    }

    #Here we fill a dataframe and a vector containing the predicted drift signals and calibration coefficient respectively for all elements
    #If i == 1, this is the first element and we need to initialize the variables
    #Else we can directly fill the dataframe and vector using dataframe[i] or vector[i]
    if (i == 1){
      driftDataFrame <- data.frame(driftPredict)
      coefVector <- summary(calibrationModel)$coefficients[1,1]
    }
    else{
      driftDataFrame[i] <- driftPredict
      coefVector[i] <- summary(calibrationModel)$coefficients[1,1]
    }
    
    #If we chose to correct for the drift (eDriftChoice[i] == TRUE)
    #Divide the drift dataframe at the ith position by the intersect at dt = 0 of the drift model
    #We transform the drift dataframe by a drift relative to dt = 0 so that we can apply the drift correction by simple multiplication
    #Else, if we chose not to correct the drift, we fill the dataframe at the ith position with 1 so that multiplication produces identity
    if (eDriftChoice[i]){
      driftDataFrame[i] <- driftDataFrame[i] / summary(driftModel)$coefficients[1,1]
      driftDataFrame[is.na(driftDataFrame[,i]), i] <- 1
    } 
    else {
      driftDataFrame[i] <- rep(1, nrow(signal[[1]]))
    }
  }
  
  signalDriftCorrectedRatio <- signal[[1]] / driftDataFrame
  signalDriftCorrectedRSD <- signal[[2]]
  concentration <- t(t(signalDriftCorrectedRatio)*coefVector)
  concentration[concentration < 0] <- "<blk"
  
  concentrationRSD <- signalDriftCorrectedRSD
  concentrationRSD[concentration < 0] <- "N/A"
  
  return(list(concentration,concentrationRSD))
}