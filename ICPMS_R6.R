library(shiny)

reactiveTrigger <- function() {
  counter <- reactiveVal(0)
  list(
    depend = function() {
      counter()
      invisible()
    },
    trigger = function() {
      counter( isolate(counter()) + 1 )
    }
  )
}

IcpDataReactive <- R6::R6Class(
  public = list(
    initialize = function(IcpData) {
      private$countValues = IcpData$countValues
      private$countSd = IcpData$countSd
    },
    getCountValues = function() {
      private$rxTrigger$depend()
      return(private$countValues)
    },
    getCountSd = function() {
      private$rxTrigger$depend()
      return(private$countSd)
    },
    getCountRsd = function() {
      private$rxTrigger$depend()
      countRsd <- private$countSd / private$countValues * 100
      return(countRsd)
    },
    getElementFullNames = function() {
      private$rxTrigger$depend()
      return( colnames(private$countValues) )
    },
    getMetadata = function() {
      private$rxTrigger$depend()
      return(private$sharedEnvironment$metadata)
    },
    getMetadataNames = function() {
      private$rxTrigger$depend()
      return( colnames(private$sharedEnvironment$metadata) ) 
    },
    setCountValuesAndSd = function(IcpData) {
      private$rxTrigger$trigger()
      private$countValues = IcpData$countValues
      private$countSd = IcpData$countSd
    },
    setCountValuesAndRsd = function(IcpData) {
      private$rxTrigger$trigger()
      private$countValues = IcpData$countValues
      private$countSd = IcpData$countRsd / 100 * IcpData$countValues
    },
    setMetadata = function(IcpMetadata) {
      private$rxTrigger$trigger()
      private$sharedEnvironment$metadata = IcpMetadata
    }
  ),
  private = list(
    sharedEnvironment = new.env(),
    countValues = NULL,
    countSd = NULL,
    rxTrigger = reactiveTrigger()
  )
)