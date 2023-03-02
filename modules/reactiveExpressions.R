#' Module responsible for computing intermediate and final results
#' 
#' @return A reactiveValues @process
#' 
#' process$ 
#' @countsPerSecond() (matrix, nrow = sampleNumber, ncol = analyteNumber + internalStandardNumber)
#' @analyteCountsPerSecond() (matrix, nrow = sampleNumber, ncol = analyteNumber)
#' @internalStandardCountsPerSecond() (matrix, nrow = sampleNumber, ncol = internalStandardNumber)
#' 
#' @relativeStandardDeviation() (matrix, nrow = sampleNumber, ncol = analyteNumber + internalStandardNumber)
#' @analyteCountsPerSecondRelativeStandardDeviation() (matrix, nrow = sampleNumber, ncol = analyteNumber)
#' @internalStandardCountsPerSecondRelativeStandardDeviation() (matrix, nrow = sampleNumber, ncol = internalStandardNumber)
#' 
#' @analyteCountsPerSecondEudc() (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' @internalStandardCountsPerSecondEudc() (eudc, nrow = sampleNumber, ncol = internalStandardNumber)
#' 
#' @internalStandardEudcAdaptedToAnalytes() (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' @analyteToIstdRatio (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' 
#' @analyteToIstdBlankRatio (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' @analyteToIstdRatioBlankCorrected (eudc, nrow = sampleNumber, ncol = analyteNumber) 
#' 
#' @driftFactorEudc (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' @driftAndBlankCorrectedEudc (eudc, nrow = sampleNumber, ncol = analyteNumber)
#' 
#' @calibrationCoefficientMatrix (double, length = analyteNumber)
#' @concentration (eudc, nrow = sampleNumber, ncol = analyteNumber)

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

reactiveExpressions_server <- function(id, fileUpload, indexCreation, interferenceCorrection, blankProcessing, driftProcessing) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      process <- reactiveValues()
      
      # Passed from fileUpload module
      extracted <- reactive(fileUpload$extracted)
      parameters <- reactive(fileUpload$parameters)
      
      # Passed from driftProcessing module
      listOfElementSpecificDriftIndex <- reactive(driftProcessing$listOfElementSpecificDriftIndex())
      elementSpecificDriftModels <- reactive(driftProcessing$elementSpecificDriftModels())
      driftCorrectedElements <- reactive(driftProcessing$driftCorrectedElements())
      
      process$countsPerSecond <- reactive({
        extracted()$main %>%
          select(ends_with("_CPS")) %>%
          as.matrix()
        # countsPerSecond <- extracted()$main[ , extracted()$secondRowOfMain == "CPS", drop=FALSE]
        # numericalCountsPerSecond <- apply(countsPerSecond, c(1,2), is.numeric)
        # countsPerSecond[!numericalCountsPerSecond] <- NA
        # as.matrix(countsPerSecond)
      })
      
      process$analyteCountsPerSecond <- reactive({
        extracted()$main %>%
          select(paste(parameters()$analyteNames, "CPS", sep = "_")) %>%
          as.matrix()
      })
      
      process$internalStandardCountsPerSecond <- reactive({
        extracted()$main %>%
          select(paste(parameters()$internalStandardNames, "CPS", sep = "_")) %>%
          as.matrix()
      })
      
      process$relativeStandardDeviation <- reactive({
        extracted()$main %>%
          select(ends_with("_CPS RSD")) %>%
          as.matrix()
        # relativeStandardDeviation <- extracted()$main[ , extracted()$secondRowOfMain == "CPS RSD", drop=FALSE]
        # numericalRelativeStandardDeviation <- apply(relativeStandardDeviation, c(1,2), is.numeric)
        # relativeStandardDeviation[!numericalRelativeStandardDeviation] <- NA
        # as.matrix(relativeStandardDeviation)
      })
      
      process$analyteCountsPerSecondRelativeStandardDeviation <- reactive({
        extracted()$main %>%
          select(paste(parameters()$analyteNames, "CPS RSD", sep = "_")) %>%
          as.matrix()
      })
      
      process$internalStandardCountsPerSecondRelativeStandardDeviation <- reactive({
        extracted()$main %>%
          select(paste(parameters()$internalStandardNames, "CPS RSD", sep = "_")) %>%
          as.matrix()
      })
      
      process$analyteCountsPerSecondEudc <- reactive({
        EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                            estimatedData = process$analyteCountsPerSecond(),
                                            uncertaintyData = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                            uncertaintyType = "rsd")
      })
      
      process$analyteCountsPerSecondBlankEudc <- reactive({
        if (input$enableBlankCorrection == TRUE)
        {
          process$analyteCountsPerSecondEudc() %>% applyModifierToEudc(blankProcessing$blankModifiers())
        }
        else
        {
          EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                              estimatedData = matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber),
                                              uncertaintyData = matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber),
                                              uncertaintyType = "sd")
        }
      })
      
      process$analyteCountsPerSecondBlankCorrected <- reactive({
        process$analyteCountsPerSecondEudc()$subtract(process$analyteCountsPerSecondBlankEudc())
      })
      
      process$interferenceCountsPerSecond <- reactive({
        process$analyteCountsPerSecondBlankCorrected() %>% applyCumulativeModifierToEudc(process$analyteCountsPerSecondEudc(), interferenceCorrection$interferenceModifiers())
      })
      
      process$interferenceCorrectedCountsPerSecond <- reactive({
        interferenceCorrected <- process$analyteCountsPerSecondEudc()$subtract(process$interferenceCountsPerSecond())
        interferenceCorrected$zeroOutNegatives()
      })
      
      process$internalStandardCountsPerSecondEudc <- reactive({
        EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$internalStandardNames,
                                            estimatedData = process$internalStandardCountsPerSecond(),
                                            uncertaintyData = process$internalStandardCountsPerSecondRelativeStandardDeviation(),
                                            uncertaintyType = "rsd")
      })
      
      process$internalStandardEudcAdaptedToAnalytes <- reactive({
        if (!is.null(extracted()$internalStandardToAnalyteAssignment)) 
        {
          return(createinternalStandardEudcAdaptedToAnalytes(internalStandardToAnalyteAssignmentDataframe = extracted()$internalStandardToAnalyteAssignment,
                                                             internalStandardCountsPerSecondEudc = process$internalStandardCountsPerSecondEudc(),
                                                             parameters = parameters()))
        }
        else
        {
          return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                                     estimatedData = matrix(1, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber),
                                                     uncertaintyData = matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber),
                                                     uncertaintyType = "sd"))
        }
      })
      
      process$analyteToIstdRatio <- reactive({
        process$interferenceCorrectedCountsPerSecond()$divideBy(process$internalStandardEudcAdaptedToAnalytes())
      })
      
      process$analyteToIstdBlankRatio <- reactive({
        if (input$enableBlankCorrection == TRUE)
        {
          process$analyteToIstdRatio() %>% applyModifierToEudc(blankProcessing$blankModifiers())
        }
        else
        {
          EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                              estimatedData = matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber),
                                              uncertaintyData = matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber),
                                              uncertaintyType = "sd")
        }
      })
      
      process$analyteToIstdRatioBlankCorrected <- reactive({
        process$analyteToIstdRatio()$subtract(process$analyteToIstdBlankRatio())
      })
      
      process$driftFactorEudc <- reactive({
        
        if (driftCorrectedElements() %>% is.null()) {
          return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                                     estimatedData = matrix(1, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber),
                                                     uncertaintyData = matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber),
                                                     uncertaintyType = "sd"))
        }
        
        for (elementIndex in seq(parameters()$analyteNumber)){
          
          elementFullName <- parameters()$analyteNames[elementIndex]
          elementSpecificDriftIndex <- listOfElementSpecificDriftIndex()[[elementFullName]]
          
          if (all(is.na(elementSpecificDriftModels()[[elementFullName]]))){
            driftPredict = rep(NA, parameters()$sampleNumber)
          }
          else{
            driftModel <- elementSpecificDriftModels()[[elementFullName]]
            driftPredict=predict(driftModel, newdata = data.frame(dt=parameters()$deltaTime()))
          }
          
          if (elementIndex == 1){
            driftMatrix <- as.matrix(driftPredict)
          }
          else{
            driftMatrix <- cbind(driftMatrix, driftPredict)
          }
          
          if (driftCorrectedElements()[elementIndex] == TRUE){
            driftMatrix[ , elementIndex] <- driftMatrix[ , elementIndex] / driftMatrix[min(elementSpecificDriftIndex), elementIndex]
            driftMatrix[is.na(driftMatrix[ , elementIndex]), elementIndex] <- 1
            driftMatrix[seq(min(elementSpecificDriftIndex)), elementIndex] <- 1
          } 
          else {
            driftMatrix[ , elementIndex] <- rep(1, parameters()$sampleNumber)
          }
        }
        return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                                   estimatedData = driftMatrix,
                                                   uncertaintyData = matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber),
                                                   uncertaintyType = "sd"))
      })
      
      process$driftAndBlankCorrectedEudc <- reactive({
        driftAndBlankCorrectedEudc <- process$analyteToIstdRatioBlankCorrected()$divideBy(process$driftFactorEudc())
        return(driftAndBlankCorrectedEudc)
      })
      
      process$calibrationCoefficientMatrix <- reactive({
        
        calibrationLinearRegressionSlope <- c()
        elementsWithCalibrationIssues <- c()
        
        for (elementIndex in 1:parameters()$analyteNumber){
          
          elementFullName <- parameters()$analyteNames[elementIndex]
          
          calibrationSignalUncertaintyConcentration <- getCalibrationData(elementFullName = elementFullName, signalMatrix = process$analyteToIstdRatioBlankCorrected(),
                                                                          standardIdentificationColumn = parameters()[["categoricalDataAndTime"]][ , "Level",drop = TRUE], standardDataFrame = extracted()$standard)
          
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
        
        return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                                   estimatedData = concentration,
                                                   uncertaintyData = concentrationRSD,
                                                   uncertaintyType = "rsd"))
      })
      
      return(list(process = process))
    }
  )
}
