# library(shiny)
# library(stringr)

agilentElementNamePattern <- "[ ]{2}[A-Z]{1}[a-z]*[ ]{2}"

createInterferenceParameters <- function(index, intereferedElement, interferingElement) {
  
  parameterClasses <- c(class(numeric)[1], class(intereferedElement)[1], class(interferingElement)[1])
  expectedClasses <- c("numeric", "character", "character")
  
  if (!identical(parameterClasses, expectedClasses)) {
    stop("The classes of the provided parameters are incorrect")
  }
  if (!identical(size(index), size(intereferedElement), size(interferingElement))) {
    stop("The character vector sizes must be the same")
  }
  
  instance <- list()
  instance$index <- index
  instance$intereferedElement <- intereferedElement
  instance$interferingElement <- interferingElement
  
  class(instance) <- "interferenceParameters"
  
  return(instance)
}

EstimationUncertaintyDataCouple <- R6Class("EstimationUncertaintyDataCouple",
  public = list(
   initialize = function(elementFullNames, estimatedData, uncertaintyData, uncertaintyType) {
     if (!is.matrix(estimatedData) | !is.matrix(uncertaintyData)) {
       stop("estimation and uncertainty data must be matrixes")
     }
     if (length(elementFullNames) != ncol(estimatedData) | length(elementFullNames) != ncol(uncertaintyData)) {
       stop("incompatible element name number and column number in estimation and/or uncertainty")
     }
     if (nrow(estimatedData) != nrow(uncertaintyData)) {
       stop("incompatible sample number between estimation and uncertainty")
     }
     private$elementFullNames <- elementFullNames
     private$estimatedData <- estimatedData
     if (uncertaintyType == "rsd") {
       private$sdData <- self$createEmptyMatrixOfEudcSize()
       notApplicableUncertaintyData <- apply(uncertaintyData, c(1,2), is.na)
       notApplicableEstimatedData <- apply(estimatedData, c(1,2), is.na)
       private$sdData[!notApplicableUncertaintyData] <- uncertaintyData[!notApplicableUncertaintyData] / 100 * estimatedData[!notApplicableUncertaintyData]
       private$sdData[notApplicableUncertaintyData & !notApplicableEstimatedData] <- 0
       private$sdData[notApplicableEstimatedData] <- NA
     }
     else if (uncertaintyType == "sd") {
       private$sdData <- uncertaintyData
     }
     colnames(private$estimatedData) <- elementFullNames
     colnames(private$sdData) <- elementFullNames
   },
   divideBy = function(otherEudc) {
     estimate <- private$estimatedData / otherEudc$getEstimation()
     uncertainty <- sqrt((self$getRsd()/100)^2 + (otherEudc$getRsd()/100)^2) * 100
     return(EstimationUncertaintyDataCouple$new(elementFullNames = private$elementFullNames, estimatedData = estimate,
            uncertaintyData = uncertainty,
            uncertaintyType = "rsd"))
   },
   multiplyBy = function(otherEudc) {
     estimate <- private$estimatedData * otherEudc$getEstimation()
     uncertainty <- sqrt((self$getRsd()/100)^2 + (otherEudc$getRsd()/100)^2) * 100
     return(EstimationUncertaintyDataCouple$new(elementFullNames = private$elementFullNames,
                                                estimatedData = estimate,
                                                uncertaintyData = uncertainty,
                                                uncertaintyType = "rsd"))
   },
   add = function(otherEudc) {
     estimate <- private$estimatedData + otherEudc$getEstimation()
     uncertainty <- sqrt(self$getSd()^2 + otherEudc$getSd()^2)
     return(EstimationUncertaintyDataCouple$new(elementFullNames = private$elementFullNames,
                                                estimatedData = estimate,
                                                uncertaintyData = uncertainty,
                                                uncertaintyType = "sd"))
   },
   subtract = function(otherEudc) {
     estimate <- private$estimatedData - otherEudc$getEstimation()
     uncertainty <- sqrt(self$getSd()^2 + otherEudc$getSd()^2)
     return(EstimationUncertaintyDataCouple$new(elementFullNames = private$elementFullNames,
                                                estimatedData = estimate,
                                                uncertaintyData = uncertainty,
                                                uncertaintyType = "sd"))
   },
   getElementFullNames = function() {
     return(private$elementFullNames)
   },
   getElementNumber = function() {
     return(length(private$elementFullNames))
   },
   getSampleNumber = function() {
     return(nrow(private$estimatedData))
   },
   getEstimation = function() {
     return(private$estimatedData)
   },
   getSd = function() {
     return(private$sdData)
   },
  getRsd = function() {
    notApplicableData <- apply(private$estimatedData, c(1,2), is.na)
    applicableEstimationEqualsToZero <- (private$estimatedData[!notApplicableData] == 0)
    rsdData <- self$createEmptyMatrixOfEudcSize()
    rsdData[!notApplicableData][!applicableEstimationEqualsToZero] <- private$sdData[!notApplicableData][!applicableEstimationEqualsToZero] / private$estimatedData[!notApplicableData][!applicableEstimationEqualsToZero] * 100
    rsdData[!notApplicableData][applicableEstimationEqualsToZero] <- NA
    rsdData[notApplicableData] <- NA
    colnames(rsdData) <- private$elementFullNames
    return(rsdData)
  },
  createEmptyMatrixOfEudcSize = function() {
    rowNumber = nrow(private$estimatedData)
    columnNumber = ncol(private$estimatedData)
    return(matrix(nrow = rowNumber, ncol = columnNumber))
  }),
   private = list(
     elementFullNames = NULL,
     estimatedData = NULL,
     sdData = NULL
   )
)

correctInterferences <- function(eudcToCorrect, interferenceParameters) {
  
  if (length(interferenceParameters$index == 0)) {
    stop("interference index of length < 1")
  }
  else if (length(interferenceParameters$index == 1)) {
    interferenceRatioValue <- eudcToCorrect$estimatedData[interferenceParameters$index, interferenceParameters$intereferedElement] / eudcToCorrect$estimatedData[interferenceParameters$index, interferenceParameters$intereferingElement]
    interferenceRatioSd <- sqrt((eudcToCorrect$rsdData[interferenceParameters$index, interferenceParameters$intereferedElement] / 100)^2 +
                                 (eudcToCorrect$rsdData[interferenceParameters$index, interferenceParameters$intereferingElement] / 100)^2) * interferenceRatioValue
  }
  else {
    interferenceRatioValues <- eudcToCorrect$estimatedData[interferenceParameters$index, interferenceParameters$intereferedElement] / eudcToCorrect$estimatedData[interferenceParameters$index, interferenceParameters$intereferingElement] %>%
    interferenceRatioValue <- mean(interferenceRatioValues)
    interferenceRatioSd <- sd(interferenceRatioValues)
  }
  interferenceValues <- interferenceRatioValue * eudcToCorrect$estimatedData[ , interferenceParameters$intereferingElement]
  interferenceSd <- sqrt((interferenceRatioSd / interferenceRatioValue)^2 + (eudcToCorrect$rsdData[ , interferenceParameters$intereferingElement] / 100)^2) * interferenceValues
  
  interferenceCorrectedValues <- eudcToCorrect$estimatedData[ , interferenceParameters$intereferedElement] - interferenceValues
  interferenceCorrectedSd <- sqrt(eudcToCorrect$sdData[ , interferenceParameters$intereferedElement]^2 + interferenceSd^2)
  
  correctedEstimatedData <- eudcToCorrect$estimatedData
  correctedEstimatedData[ , interferenceParameters$intereferedElement] <- interferenceCorrectedValues
  
  correctedEstimatedData <- eudcToCorrect$estimatedData
  
  correctedUncertaintyData <- 
  correctedUncertaintyType <- "sd"
  correctedEudc <- createEstimationUncertaintyDataCouple()
  
  return(correctedEudc)
}

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
    eSignal <- signalMatrix$getEstimation()[standardRowNumberInMain, elementFullName]
    eRSD <- signalMatrix$getSd()[standardRowNumberInMain, elementFullName]
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

replaceValues <- function(eudcToModify, howToReplace, linesToBeReplaced, linesUsedForReplacement, parameters){
  
  signalEstimation = eudcToModify$getEstimation()
  signalRelativeStandardDeviation = eudcToModify$getRsd()
  elementFullNames <- eudcToModify$getElementFullNames()
  modifiedEstimation <- signalEstimation
  modifiedRelativeStandardDeviation <- signalRelativeStandardDeviation
  
  #This condition is useful in shiny tables to return the table signal and RSD when no replacements are required
  if (howToReplace == "none" | is.null(linesToBeReplaced)) {
    return(eudcToModify)
  }
  
  lineNumber <- eudcToModify$getSampleNumber()
  
  # Forces logical index
  if (!is.logical(linesToBeReplaced)) {
    linesToBeReplaced <- convertToLogical(linesToBeReplaced, lineNumber)
  }
  # The function accepts a NULL sourceLineIndex. The sourceLineIndex variable takes the opposite value of the logical lineIndex (!lineIndex)
  if (is.null(linesUsedForReplacement)){
    linesUsedForReplacement <- !linesToBeReplaced
  }
  # Forces logical index
  else if (!is.logical(linesUsedForReplacement)) {
    linesUsedForReplacement <- convertToLogical(linesUsedForReplacement, lineNumber)
  }
  # If linesUsedForReplacement is the "all" index, converts linesUsedForReplacement to the inverse of the linesToBeReplaced
  if (identical(linesUsedForReplacement, rep(TRUE, rep = lineNumber))) {
    linesUsedForReplacement <- !linesToBeReplaced
  }

  numericalLinesToBeReplaced <- which(linesToBeReplaced)
  numericalLinesUsedForReplacement <- which(linesUsedForReplacement)
  
  for (elementFullName in elementFullNames){
    
    for (lineToReplace in numericalLinesToBeReplaced){
      
      previousLine <- getClosestInIndex(lineToReplace, numericalLinesUsedForReplacement, "prev")
      previousEstimation <- signalEstimation[previousLine, elementFullName]
      previousRelativeStandardDeviation <- signalRelativeStandardDeviation[previousLine, elementFullName]

      if (howToReplace == "mean"){
        
        nextLine <- getClosestInIndex(lineToReplace, numericalLinesUsedForReplacement, "next")
        nextEstimation <- signalEstimation[nextLine, elementFullName]
        
        if (is.na(previousEstimation) | is.na(nextEstimation) | is.na(modifiedEstimation[lineToReplace, elementFullName])) {
          modifiedEstimation[lineToReplace, elementFullName] <- NA
        }
        else {
          modifiedEstimation[lineToReplace, elementFullName] <- mean(c(previousEstimation, nextEstimation))
          modifiedRelativeStandardDeviation[lineToReplace,elementFullName] <- calculateRSD(c(previousEstimation, nextEstimation))
        }
        
      } 
      else if (howToReplace == "prev"){
        
        if (is.na(previousEstimation) | is.na(modifiedEstimation[lineToReplace,elementFullName])) {
          modifiedEstimation[lineToReplace, elementFullName] <- NA
        }
        else {
          modifiedEstimation[lineToReplace, elementFullName] <- previousEstimation
          modifiedRelativeStandardDeviation[lineToReplace, elementFullName] <- previousRelativeStandardDeviation
        }
        
      }
      else if (howToReplace == "averageInBlankIndex"){
        
        if (length(numericalLinesUsedForReplacement) > 1) {
          modifiedEstimation[lineToReplace, elementFullName] <- mean(signalEstimation[numericalLinesUsedForReplacement, elementFullName])
          modifiedRelativeStandardDeviation[lineToReplace,elementFullName] <- sd(signalEstimation[numericalLinesUsedForReplacement, elementFullName]) / signalEstimation[lineToReplace, elementFullName] * 100
        }
        else if (length(numericalLinesUsedForReplacement) == 1) 
        {
          modifiedEstimation[lineToReplace,elementFullName] <- signalEstimation[numericalLinesUsedForReplacement, elementFullName]
          modifiedRelativeStandardDeviation[lineToReplace,elementFullName] <- signalRelativeStandardDeviation[numericalLinesUsedForReplacement, elementFullName]
        }
        else 
        {
          modifiedEstimation[lineToReplace, elementFullName] <- NA
          modifiedRelativeStandardDeviation[lineToReplace, elementFullName] <- NA
        }
        
      }
    }
  }
  return(EstimationUncertaintyDataCouple$new(elementFullNames = elementFullNames,
                                             estimatedData = modifiedEstimation,
                                             uncertaintyData = modifiedRelativeStandardDeviation,
                                             uncertaintyType = "rsd"))
}

##Function to create the ISTD template based on raw datafile
##Returns a dataframe containing a column analyte with the analyte names
##and a column ISTD with the ISTD names
createISTDtemplate <- function(dataFileName, sep){
  
  #header_1 contains the first line of the file
  header_1 <- scan(dataFileName, nlines = 1, what = character(),sep = sep)
  #header_2 contains the second line of the file
  header_2 <- scan(dataFileName, nlines = 1, what = character(),sep = sep, skip=1)
  
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
  
  
createinternalStandardEudcAdaptedToAnalytes <- function(internalStandardToAnalyteAssignmentDataframe, internalStandardCountsPerSecondEudc, parameters){
  
  
  analyteNames <- parameters$analyteNames
  istdColumnPosition <- 2
  
  estimationMatrix <- matrix(nrow = parameters$sampleNumber, ncol = parameters$analyteNumber)
  uncertaintyMatrix <- matrix(nrow = parameters$sampleNumber, ncol = parameters$analyteNumber)
  
  for (analyteIncrement in seq(1, parameters$analyteNumber)){
    analyteSpecificIstd <- internalStandardToAnalyteAssignmentDataframe[analyteIncrement, istdColumnPosition]
    estimationMatrix[ , analyteIncrement] <- internalStandardCountsPerSecondEudc$getEstimation()[ , analyteSpecificIstd]
    uncertaintyMatrix[ , analyteIncrement] <- internalStandardCountsPerSecondEudc$getSd()[, analyteSpecificIstd]
  }
  
  colnames(estimationMatrix) <- analyteNames
  colnames(uncertaintyMatrix) <- analyteNames
  
  return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                             estimatedData = estimationMatrix,
                                             uncertaintyData = uncertaintyMatrix,
                                             uncertaintyType = "sd"))
}