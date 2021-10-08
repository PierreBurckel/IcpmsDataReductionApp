# library(shiny)
# library(stringr)

agilentElementNamePattern <- "[ ]{2}[A-Z]{1}[a-z]*[ ]{2}"

# setMatrixFormat <- function(matrixToFormat, matrixToFormatColumnNames = NULL, columnsToAdd, columnsToAddColumnNames = NULL) {
#   
#   sampleNumber <- nrow(matrixToFormat)
#   
#   if (class(matrixToFormat) == "matrix" || class(matrixToFormat) ==  "data.frame") {
#     matrixToFormat <- as.matrix(matrixToFormat)
#     formatedMatrix <- cbind(columnsToAdd, matrixToFormat)
#     colnames(formatedMatrix) <- c(columnsToAddColumnNames, matrixToFormatColumnNames)
#   }
#   
#   if (class(matrixToFormat) == "list") {
#     matrixToFormat[[1]] <- as.matrix(matrixToFormat[[1]])
#     matrixToFormat[[2]] <- as.matrix(matrixToFormat[[2]])
#     formatedMatrix <- matrix(0, nrow = sampleNumber, ncol = 2 * analyteNumber)
#     formatedMatrixColumnNames <- rep("", 2 * analyteNumber)
#     
#     firstMatrixInsertionIndex <- seq(from = 1, to = 2 * analyteNumber, by = 2)
#     secondMatrixInsertionIndex <- seq(from = 2, to = 2 * analyteNumber, by = 2)
#     
#     formatedMatrix[ , firstMatrixInsertionIndex] <- matrixToFormat[[1]]
#     formatedMatrix[ , secondMatrixInsertionIndex] <- matrixToFormat[[2]]
#     formatedMatrixColumnNames[firstMatrixInsertionIndex] <- matrixToFormatColumnNames
#     formatedMatrixColumnNames[secondMatrixInsertionIndex] <- paste0(matrixToFormatColumnNames, " RSD (%)")
#     
#     formatedMatrix <- cbind(columnsToAdd, formatedMatrix)
#     colnames(formatedMatrix) <- c(columnNamesToAdd, formatedMatrixColumnNames)
#   }
#   
#   return(formatedMatrix)
# }

createStandardDataFrameFromFile <- function(dataPath, sep) {
  
  standardDataFrame <- read.table(dataPath, header = FALSE, sep = sep, stringsAsFactors=FALSE)
  
  if (!identical(unique(standardDataFrame[1, ]), standardDataFrame[1, ])) {
    shinyalert("Impossible to extract", "In the standard file, the standard labels in the first row are not unique", type = "error")
    return(NULL)
  }
  
  if (!identical(unique(standardDataFrame[ ,1]), standardDataFrame[ ,1])) {
    shinyalert("Impossible to extract", "In the standard file, the analyte names in the first column are not unique", type = "error")
    return(NULL)
  }
  
  colnames(standardDataFrame) <- standardDataFrame[1, ]
  standardDataFrame <- standardDataFrame[-1, ]
  row.names(standardDataFrame) <- standardDataFrame[ , 1]
  standardDataFrame <- standardDataFrame[ , -1]
  
  return(standardDataFrame)
}

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

isAnalyteValidInStandardDataFrame <- function(elementFullName, standardDataFrame) {
  
  if (!(elementFullName %in% row.names(standardDataFrame))) {
    shinyalert("Issue with standard file",
               paste0(elementFullName, " not found in the standard file. If not in the standards, create an empty line for ", elementFullName, " in the file."),
               type = "error")
    return(FALSE)
  }
  
  elementSpecificStandardIndex <- which(!is.na(standardDataFrame[elementFullName, ]))
  
  if (length(elementSpecificStandardIndex) == 0) {
    shinyalert("Issue with standard file",
               paste0("The standard file line for ",elementFullName, " should be of size > 0"),
               type = "error")
    return(FALSE)
  }
  
  return(TRUE)
}

getElementSpecificDriftIndex <- function(elementFullName, standardDataFrame, standardIdentificationColumn, driftIndex) {
  
  if (isAnalyteValidInStandardDataFrame(elementFullName, standardDataFrame) == FALSE) {
    return(NULL)
  }
  
  elementSpecificStandardIndex <- which(!is.na(standardDataFrame[elementFullName, ]))
  
  standardIdColumnInData <- make.names(standardIdentificationColumn)
  elementSpecificStandardIdHeaderInStandardFile <- make.names(colnames(standardDataFrame)[elementSpecificStandardIndex])
  isStandardBothInDataAndStandardFile <- standardIdColumnInData %in% elementSpecificStandardIdHeaderInStandardFile
  firstStandardNumIndex <- min(which(isStandardBothInDataAndStandardFile))
  elementSpecificDritftIndex <- driftIndex[driftIndex >= firstStandardNumIndex]
    
  return(elementSpecificDritftIndex)
}

getCalibrationData <- function(elementFullName, signalMatrix, standardIdentificationColumn, standardDataFrame) {
  
  if (!(elementFullName %in% row.names(standardDataFrame))) {
    return(NULL)
  }
  
  standardNamesInStandardFile <- colnames(standardDataFrame[elementFullName, ])
  standardConcentrationInStandardFile <- standardDataFrame[elementFullName, ]
  
  standardRowNumberInMain <- as.numeric(sapply(paste("^", standardNamesInStandardFile, "$", sep=""),
                                   grep, standardIdentificationColumn))
  
  isInMainAndStandardFiles <- !is.na(standardRowNumberInMain) & !is.na(standardConcentrationInStandardFile)
  
  standardConcentrationInStandardFile <- standardConcentrationInStandardFile[isInMainAndStandardFiles]
  standardRowNumberInMain <- standardRowNumberInMain[isInMainAndStandardFiles]
  
  if (length(standardRowNumberInMain) == 0 || length(standardConcentrationInStandardFile) == 0)
  {
    return(NULL)
  }
  else 
  {
    eSignal <- signalMatrix[[1]][standardRowNumberInMain, elementFullName]
    eRSD <- signalMatrix[[2]][standardRowNumberInMain, elementFullName]
  }
  
  return(cbind(Signal=eSignal, RSD=eRSD, Concentration=as.numeric(standardConcentrationInStandardFile)))
}

getElementName <- function(elementFullName, pattern) {
  return(gsub(" ", "", str_extract(elementFullName, pattern), fixed = TRUE))
}

mergeMatrixes <- function(matrix1, matrix2) {
  
  if (!identical(dim(matrix1), dim(matrix2)))
  {
    stop("Merged matrixes need to have the same dimensions")
  }
  
  rowNumber <- nrow(matrix1)
  columnNumber <- ncol(matrix1)
  
  matrix1 <- as.matrix(matrix1)
  matrix2 <- as.matrix(matrix2)
  mergedMatrix <- matrix(0, nrow = rowNumber, ncol = 2 * columnNumber)
  
  firstMatrixInsertionIndex <- seq(from = 1, to = 2 * columnNumber, by = 2)
  secondMatrixInsertionIndex <- seq(from = 2, to = 2 * columnNumber, by = 2)
  
  mergedMatrix[ , firstMatrixInsertionIndex] <- matrix1
  mergedMatrix[ , secondMatrixInsertionIndex] <- matrix2
  
  return(mergedMatrix)
}

createColoredTextForBooleanViewing <- function(state, stateText, invalidStateText, validStateText) {
  
  if (state == TRUE) {
    displayedText <- paste('<span style=\"color:green\">', stateText, validStateText)
  }
  else if (state == FALSE) {
    displayedText <- paste('<span style=\"color:red\">', stateText, invalidStateText)
  }
  return(displayedText)
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
    if (mean(values) == 0) {
      return(0)
    } else {
      return(sd(values)/mean(values)*100)
    }
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
replaceValues <- function(sourceValues, colRange, replaceType, lineIndex, sourceLineIndex){
  
  ICPsignal = sourceValues[["signal"]]
  RSD = sourceValues[["RSD"]]
  modifiedIcpSignal <- ICPsignal
  modifiedIcpRsd <- RSD
  
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
        modifiedIcpSignal[nLine,j] <- mean(c(prevValue,nextValue))
        modifiedIcpRsd[nLine,j] <- calculateRSD(c(prevValue,nextValue))
      } 
      else if (replaceType == "prev"){
        modifiedIcpSignal[nLine,j] <- prevValue
        modifiedIcpRsd[nLine,j] <- prevValueRSD
      }
      else if (replaceType == "averageInBlankIndex"){
        if (length(nSourceLineIndex) > 1) {
          modifiedIcpSignal[nLine,j] <- mean(ICPsignal[nSourceLineIndex, j])
          modifiedIcpRsd[nLine,j] <- sd(ICPsignal[nSourceLineIndex, j]) / ICPsignal[nLine,j] * 100
        }
        else if (length(nSourceLineIndex) == 1) 
        {
          modifiedIcpSignal[nLine,j] <- ICPsignal[nSourceLineIndex, j]
          modifiedIcpRsd[nLine,j] <- RSD[nSourceLineIndex, j]
        }
        else 
        {
          modifiedIcpSignal[nLine,j] <- NA
          modifiedIcpRsd[nLine,j] <- NA
          
        }
      }
    }
  }
  return(list(signal=modifiedIcpSignal,RSD=modifiedIcpRsd))
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
  standardDataFrame <- read.table(stdFileName, header = TRUE, sep=';', stringsAsFactors=FALSE)
  
  #Remove potential duplicates of element names in the standard dataframe
  standardDataFrame <- removeDuplicateLines(standardDataFrame)
  
  #Assigns the first column of the std file to the row names then removes the first column
  row.names(standardDataFrame) <- standardDataFrame[,1]
  standardDataFrame <- standardDataFrame[,2:length(standardDataFrame)]
  
  return(list(raw=dat.raw, std=standardDataFrame, header_1=header_1, header_2=header_2))
}
  
  
createInternalStandardMatrixAdaptedToAnalytes <- function(ISTD_file, ISTD){
  
  if(is.null(ISTD_file) | is.null(ISTD)){return(1)}
  
  ISTD.matrix <- list(signal = matrix(nrow = nrow(ISTD[["signal"]]), ncol = nrow(ISTD_file)),
                      RSD = matrix(nrow = nrow(ISTD[["signal"]]), ncol = nrow(ISTD_file)))
  
  for (j in 1:nrow(ISTD_file)){
    eISTD <- ISTD_file[j,2]
    ISTD.matrix[["signal"]][ , j] <- ISTD[["signal"]][, eISTD]
    ISTD.matrix[["RSD"]][ , j] <- ISTD[["RSD"]][, eISTD]
  }
  
  colnames(ISTD.matrix[["signal"]]) <- ISTD_file[ ,1]
  colnames(ISTD.matrix[["RSD"]]) <- ISTD_file[ ,1]
  
  return(ISTD.matrix)
}