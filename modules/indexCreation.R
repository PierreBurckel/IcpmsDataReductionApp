#' Module responsible for creating line indexes
#' 
#' @return A reactiveValues @rowIndexInMain
#' 
#' rowIndexInMain$ @custom (list, int) @index_rowsMatchingRegularExpression (int)

indexCreation_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Index creation",
           sidebarLayout(
             sidebarPanel(
               selectInput(ns("searchIndexwhere"), "In:", c("Sample Names" = "smp", "Levels" = "lvl", "Type" = "type")),
               textInput(ns("searchIndexwhat"), "Search:", ""),
               selectInput(ns("searchIndexhow"), "Using:", c("Regular expression" = "regexp", "Exact match" = "ematch")),
               textInput(ns("searchIndexName"), "Index Name:", ""),
               actionButton(ns("indexSelectAll"), "Select All"),
               actionButton(ns("searchIndexCreate"), "Create new index")),
             mainPanel(
               selectInput(ns("searchIndexDisplay"), "Show:", c("Internal standards" = "ISTD", "Analytes" = "analytes")),
               DT::DTOutput(ns("indexTable"))
             )
           )
  ) 
}

indexCreation_server <- function(id, fileUpload, reactiveExpressions) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      rowIndexInMain <- reactiveValues()
      
      parameters <- reactive({fileUpload$parameters})
      process <- reactive({reactiveExpressions$process})
      
      indexTable <- reactive({
        searchWhere = input$searchIndexwhere
        firstColumnName = ""
        if (searchWhere == 'smp'){
          headerRows = parameters()[["categoricalDataAndTime"]][ , "Sample Name", drop = TRUE]
          firstColumnName = "Sample Name"}
        else if (searchWhere == 'lvl'){
          headerRows = parameters()[["categoricalDataAndTime"]][ , "Level", drop = TRUE]
          firstColumnName = "Level"}
        else if (searchWhere == 'type'){
          headerRows = parameters()[["categoricalDataAndTime"]][ , "Type", drop = TRUE]
          firstColumnName = "Type"}
        searchWhat = input$searchIndexwhat
        searchType = input$searchIndexhow
        displayWhat = input$searchIndexDisplay
        indexTable <- cbind(headerRows, process()$internalStandardCountsPerSecondEudc()$getEstimation())
        if (searchType == 'ematch'){
          searchWhat = paste("^", searchWhat, "$", sep="")
        }
        if (searchWhat == ""){
          rowIndexInMain$index_rowsMatchingRegularExpression <- seq(parameters()$sampleNumber)
        }
        else {
          rowIndexInMain$index_rowsMatchingRegularExpression <- tryCatch(
            {
              grep(searchWhat, headerRows)
            },
            error=function(cond) {
              message("Error, impossible to evaluate regular expression")
              return(NULL)
            },
            warning=function(cond) {
              message("Warning raised evaluating regular expression")
              return(NULL)
            }
          )
          # rowIndexInMain$index_rowsMatchingRegularExpression <- grep(searchWhat, headerRows)
        }
        if (displayWhat == "ISTD" & !is.null(process()$internalStandardCountsPerSecondEudc()$getEstimation())){
          indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process()$internalStandardCountsPerSecondEudc()$getEstimation()[rowIndexInMain$index_rowsMatchingRegularExpression,,drop = FALSE])
          colnames(indexTable) <- c(firstColumnName, parameters()$internalStandardNames)
        }
        else if (displayWhat == "analytes" & !is.null(process()$analyteCountsPerSecondEudc()$getEstimation())){
          indexTable = cbind(headerRows[rowIndexInMain$index_rowsMatchingRegularExpression], process()$analyteCountsPerSecondEudc()$getEstimation()[rowIndexInMain$index_rowsMatchingRegularExpression,,drop = FALSE])
          colnames(indexTable) <- c(firstColumnName, parameters()$analyteNames)
        }
        indexTable
      })
      # indexTable <- tibble()
      
      observeEvent(parameters()$sampleNumber, {
        rowIndexInMain$custom[["All"]] <- which(rep(x = TRUE, parameters()$sampleNumber))
      })
      
      output$indexTable <- DT::renderDT(datatable({
        indexTable()
      }, extensions = c('Scroller', 'Buttons'),
      options = list(dom = 'Bt', ordering=F, autoWidth = TRUE,
                     scrollX = TRUE, scrollY = 300, deferRender = TRUE, scroller = TRUE,
                     buttons = c('copy', 'csv'),
                     columnDefs = list(list(width = '120px', targets = "_all")))))
      
      indexTab_selectedRows <- reactive({
        sort(input$indexTable_rows_selected)
      })
      
      observeEvent(input$searchIndexwhat, {
        updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
      })
      
      observeEvent(input$searchIndexCreate, {
        indexName = input$searchIndexName
        customNumIndex = c(1:parameters()$sampleNumber)[rowIndexInMain$index_rowsMatchingRegularExpression][indexTab_selectedRows()]
        if (is.null(rowIndexInMain$custom)) {
          rowIndexInMain$custom <- list()
        }
        rowIndexInMain$custom[[input$searchIndexName]] <- customNumIndex
      })
      
      indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
      
      observeEvent(input$indexSelectAll, {
        if (identical(indexTab_selectedRows(), seq(nrow(indexTable()))))
        {
          DT::selectRows(indexTableProxy, NULL)
        }
        else
        {
          DT::selectRows(indexTableProxy, seq(nrow(indexTable())))
        }
        
      })
      
      observeEvent(rowIndexInMain$custom, {
        if (is.null(rowIndexInMain$custom)){return()}
        updateSelectInput(session,"sliderInput_InterferenceTab_indexForCorrectionFactorCalculation", label  = "Compute correction factor from:", choices=names(rowIndexInMain$custom),names(rowIndexInMain$custom)[1])
        updateSelectInput(session,"viewConcentrationIndex", label  = "View index:", choices=c("All", names(rowIndexInMain$custom)),"All")
      })
      
      return(rowIndexInMain)
    }
  )
}