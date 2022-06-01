reactiveExpressions_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Processing parameters",
           fluidPage(
             fluidRow(
               column(width = 4,
                      checkboxInput(ns("enableBlankCorrection"), "Enable blank correction", value = TRUE))
             )
           )
  )
}

reactiveExpressions_server <- function(id, fileUpload, indexCreation, blankProcessing) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      process <- reactiveValues()

      extracted = fileUpload$extracted
      parameters = fileUpload$parameters
      
      process$countsPerSecond <- reactive({
        countsPerSecond <- extracted$main[ , extracted$secondRowOfMain == "CPS", drop=FALSE]
        numericalCountsPerSecond <- apply(countsPerSecond, c(1,2), is.numeric)
        countsPerSecond[!numericalCountsPerSecond] <- NA
        as.matrix(countsPerSecond)
      })
      
      process$analyteCountsPerSecond <- reactive({
        process$countsPerSecond()[ , parameters$analyteNames, drop=FALSE]
      })
      
      process$internalStandardCountsPerSecond <- reactive({
        process$countsPerSecond()[ , parameters$internalStandardNames, drop=FALSE]
      })
      
      process$relativeStandardDeviation <- reactive({
        relativeStandardDeviation <- extracted$main[ , extracted$secondRowOfMain == "CPS RSD", drop=FALSE]
        numericalRelativeStandardDeviation <- apply(relativeStandardDeviation, c(1,2), is.numeric)
        relativeStandardDeviation[!numericalRelativeStandardDeviation] <- NA
        as.matrix(relativeStandardDeviation)
      })
      
      process$analyteCountsPerSecondRelativeStandardDeviation <- reactive({
        process$relativeStandardDeviation()[ , parameters$analyteNames, drop=FALSE]
      })
      
      process$internalStandardCountsPerSecondRelativeStandardDeviation <- reactive({
        process$relativeStandardDeviation()[ , parameters$internalStandardNames, drop=FALSE]
      })
      
      process$analyteCountsPerSecondEudc <- reactive({
        EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                            estimatedData = process$analyteCountsPerSecond(),
                                            uncertaintyData = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                            uncertaintyType = "rsd")
      })
      
      process$internalStandardCountsPerSecondEudc <- reactive({
        EstimationUncertaintyDataCouple$new(elementFullNames = parameters$internalStandardNames,
                                            estimatedData = process$internalStandardCountsPerSecond(),
                                            uncertaintyData = process$internalStandardCountsPerSecondRelativeStandardDeviation(),
                                            uncertaintyType = "rsd")
      })
      
      process$internalStandardEudcAdaptedToAnalytes <- reactive({
        if (!is.null(extracted$internalStandardToAnalyteAssignment)) 
        {
          return(createinternalStandardEudcAdaptedToAnalytes(internalStandardToAnalyteAssignmentDataframe = extracted$internalStandardToAnalyteAssignment,
                                                             internalStandardCountsPerSecondEudc = process$internalStandardCountsPerSecondEudc(),
                                                             parameters = parameters))
        }
        else
        {
          return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                     estimatedData = matrix(1, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                                     uncertaintyData = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                                     uncertaintyType = "sd"))
        }
      })
      
      process$analyteToIstdRatio <- reactive({
        process$analyteCountsPerSecondEudc()$divideBy(process$internalStandardEudcAdaptedToAnalytes())
      })
      
      process$analyteToIstdBlankRatio <- reactive({
        if (input$enableBlankCorrection == TRUE)
        {
          process$analyteToIstdRatio() %>% applyModifierToEudc(blankProcessing$blankModifiers())
        }
        else
        {
          EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                              estimatedData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                              uncertaintyData = matrix(0, ncol = parameters$analyteNumber, nrow = parameters$sampleNumber),
                                              uncertaintyType = "sd")
        }
      })
      
      process$analyteToIstdRatioBlankCorrected <- reactive({
        process$analyteToIstdRatio()$subtract(process$analyteToIstdBlankRatio())
      })
      
      process$elementSpecificDriftModels <- reactive({
        
        listOfDriftModels <- list()
        
        for (elementIndex in seq(parameters$analyteNumber)){
          
          elementFullName <- parameters$analyteNames[elementIndex]
          
          elementSpecificDriftIndex <- parameters$listOfElementSpecificDriftIndex()[[elementFullName]]
          
          driftSignalAndDeltaTime <- cbind(Signal=process$analyteToIstdRatioBlankCorrected()$getEstimation()[elementSpecificDriftIndex,elementIndex],
                                           dt=parameters$deltaTime()[elementSpecificDriftIndex])  
          
          if (length(elementSpecificDriftIndex) < 3){
            
            driftModel = NA
          }
          else{
            driftModel <- lm(Signal ~ poly(dt, degree=2, raw=TRUE), data = as.data.frame(driftSignalAndDeltaTime))
          }
          
          listOfDriftModels[[elementFullName]] <- driftModel
        }
        
        return(listOfDriftModels)
      })
      
      process$driftFactorEudc <- reactive({
        for (elementIndex in seq(parameters$analyteNumber)){
          
          elementFullName <- parameters$analyteNames[elementIndex]
          elementSpecificDriftIndex <- parameters$listOfElementSpecificDriftIndex()[[elementFullName]]
          
          if (all(is.na(process$elementSpecificDriftModels()[[elementFullName]]))){
            driftPredict = rep(NA, parameters$sampleNumber)
          }
          else{
            driftModel <- process$elementSpecificDriftModels()[[elementFullName]]
            driftPredict=predict(driftModel, newdata = data.frame(dt=parameters$deltaTime()))
          }
          
          if (elementIndex == 1){
            driftMatrix <- as.matrix(driftPredict)
          }
          else{
            driftMatrix <- cbind(driftMatrix, driftPredict)
          }
          
          if (parameters$driftCorrectedElements[elementIndex] == TRUE){
            driftMatrix[ , elementIndex] <- driftMatrix[ , elementIndex] / driftMatrix[min(elementSpecificDriftIndex), elementIndex]
            driftMatrix[is.na(driftMatrix[ , elementIndex]), elementIndex] <- 1
            driftMatrix[seq(min(elementSpecificDriftIndex)), elementIndex] <- 1
          } 
          else {
            driftMatrix[ , elementIndex] <- rep(1, parameters$sampleNumber)
          }
        }
        return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                   estimatedData = driftMatrix,
                                                   uncertaintyData = matrix(0, nrow = parameters$sampleNumber, ncol = parameters$analyteNumber),
                                                   uncertaintyType = "sd"))
      })
      
      process$driftAndBlankCorrectedEudc <- reactive({
        driftAndBlankCorrectedEudc <- process$analyteToIstdRatioBlankCorrected()$divideBy(process$driftFactorEudc())
        return(driftAndBlankCorrectedEudc)
      })
      
      process$calibrationCoefficientMatrix <- reactive({
        
        calibrationLinearRegressionSlope <- c()
        elementsWithCalibrationIssues <- c()
        
        for (elementIndex in 1:parameters$analyteNumber){
          
          elementFullName <- parameters$analyteNames[elementIndex]
          
          calibrationSignalUncertaintyConcentration <- getCalibrationData(elementFullName = elementFullName, signalMatrix = process$analyteToIstdRatioBlankCorrected(),
                                                                          standardIdentificationColumn = parameters[["categoricalDataAndTime"]][ , "Level"], standardDataFrame = extracted$standard)
          
          if (!is.null(calibrationSignalUncertaintyConcentration)) {
            calibrationModel <- lm(Concentration ~ 0+Signal, data=as.data.frame(calibrationSignalUncertaintyConcentration))
          }
          
          if (is.null(calibrationSignalUncertaintyConcentration) || dim(summary(calibrationModel)$coefficients)[1] == 0) {
            calibrationLinearRegressionSlope <- c(calibrationLinearRegressionSlope, NA)
            elementsWithCalibrationIssues <- c(elementsWithCalibrationIssues, elementFullName)
          }
          else
          {
            calibrationLinearRegressionSlope <- c(calibrationLinearRegressionSlope, summary(calibrationModel)$coefficients[1,1])
          }
          
        }
        
        if (!is.null(elementsWithCalibrationIssues)) {
          shinyalert("Calibration alert", paste("Impossible to compute calibration for", paste(elementsWithCalibrationIssues, collapse = ' '), sep = " "), type = "error")
        }
        
        return(calibrationLinearRegressionSlope)
      })
      
      process$concentration <- reactive({
        
        concentration <- t(t(process$driftAndBlankCorrectedEudc()$getEstimation())*process$calibrationCoefficientMatrix())
        concentrationRSD <- process$driftAndBlankCorrectedEudc()$getRsd()
        
        return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters$analyteNames,
                                                   estimatedData = concentration,
                                                   uncertaintyData = concentrationRSD,
                                                   uncertaintyType = "rsd"))
      })
      
      return(list(process = process))
    }
  )
}
