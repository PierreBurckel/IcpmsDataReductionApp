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

reactiveExpressions_server <- function(id, fileUpload, indexCreation, blankProcessing, driftProcessing) {
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
        countsPerSecond <- extracted()$main %>% select(ends_with("_CPS"))
        # countsPerSecond <- extracted()$main[ , extracted()$secondRowOfMain == "CPS", drop=FALSE]
        # numericalCountsPerSecond <- apply(countsPerSecond, c(1,2), is.numeric)
        # countsPerSecond[!numericalCountsPerSecond] <- NA
        # as.matrix(countsPerSecond)
      })
      
      process$analyteCountsPerSecond <- reactive({
        browser()
        process$countsPerSecond() %>% select(paste(parameters()$analyteNames, "CPS", sep = "_"))
      })
      
      process$internalStandardCountsPerSecond <- reactive({
        process$countsPerSecond() %>% select(paste(parameters()$internalStandardNames, "CPS", sep = "_"))
      })
      
      process$relativeStandardDeviation <- reactive({
        relativeStandardDeviation <- extracted()$main %>% select(ends_with("_CPS RSD"))
        # relativeStandardDeviation <- extracted()$main[ , extracted()$secondRowOfMain == "CPS RSD", drop=FALSE]
        # numericalRelativeStandardDeviation <- apply(relativeStandardDeviation, c(1,2), is.numeric)
        # relativeStandardDeviation[!numericalRelativeStandardDeviation] <- NA
        # as.matrix(relativeStandardDeviation)
      })
      
      process$analyteCountsPerSecondRelativeStandardDeviation <- reactive({
        process$relativeStandardDeviation() %>% select(paste(parameters()$analyteNames, "CPS RSD", sep = "_"))
      })
      
      process$internalStandardCountsPerSecondRelativeStandardDeviation <- reactive({
        process$relativeStandardDeviation() %>% select(paste(parameters()$internalStandardNames, "CPS RSD", sep = "_"))
      })
      
      process$analyteCountsPerSecondEudc <- reactive({
        EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                            estimatedData = process$analyteCountsPerSecond(),
                                            uncertaintyData = process$analyteCountsPerSecondRelativeStandardDeviation(),
                                            uncertaintyType = "rsd")
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
                                                     estimatedData = as_tibble(matrix(1, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber)),
                                                     uncertaintyData = as_tibble(matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber)),
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
          EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                              estimatedData = as_tibble(matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber)),
                                              uncertaintyData = as_tibble(matrix(0, ncol = parameters()$analyteNumber, nrow = parameters()$sampleNumber)),
                                              uncertaintyType = "sd")
        }
      })
      
      process$analyteToIstdRatioBlankCorrected <- reactive({
        process$analyteToIstdRatio()$subtract(process$analyteToIstdBlankRatio())
      })
      
      process$driftFactorEudc <- reactive({
        
        if (driftCorrectedElements() %>% is.null()) {
          return(EstimationUncertaintyDataCouple$new(elementFullNames = parameters()$analyteNames,
                                                     estimatedData = as_tibble(matrix(1, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber)),
                                                     uncertaintyData = as_tibble(matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber)),
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
                                                   estimatedData = as_tibble(driftMatrix),
                                                   uncertaintyData = as_tibble(matrix(0, nrow = parameters()$sampleNumber, ncol = parameters()$analyteNumber)),
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
                                                                          standardIdentificationColumn = parameters()[["categoricalDataAndTime"]][ , "Level"], standardDataFrame = extracted()$standard)
          
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
