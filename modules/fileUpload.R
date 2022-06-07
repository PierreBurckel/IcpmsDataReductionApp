fileUpload_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "File upload and parameters",
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("uploadedFile"), "File upload"),
        actionButton(ns("setAsMainFile"), "Set as main file"),
        actionButton(ns("setAsStandardFile"), "Set as standard file"),
        downloadButton(ns("downloadISTDTemplate"), "Download ISTD template"),
        actionButton(ns("setAsInternalStandardFile"), "Set as ISTD file"),
        actionButton(ns("extract"), "Extract"),
        htmlOutput(ns("extractionReadyText")),
        selectInput(inputId = ns("csvDelimitation"), label = "Select csv field separator (the separator will be used during extraction,
                                 you therefore need to use the same separator for all your files)",
                    choices = c(";", ","), selected = ";"),
        textInput(inputId = ns("dateFormat"), label = "Date format", value = "%d/%m/%Y %H:%M")),
      mainPanel(
        htmlOutput(ns("mainFileAssignmentText")),
        htmlOutput(ns("standardFileAssignmentText")),
        htmlOutput(ns("internalStandardFileAssignmentText")),
        div(style="display:inline-block", numericInput(inputId="fileUpload_nrow", label="Rows", value = 10, min = 1)),
        div(style="display:inline-block", numericInput(inputId="fileUpload_ncolumn", label="Columns", value = 6, min = 1)),
        tableOutput(ns("uploadedFilePreviewTable"))
      )
    )
  )
}

fileUpload_server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      uploadedFile <- reactiveValues()
      
      extracted <- reactiveValues()
      parameters <- reactiveValues()
      applicationState <- reactiveValues(isExtractionSuccessful = FALSE)
      
      parameters$deltaTime <- reactive({
        deltaTime <- as.numeric(parameters[["categoricalDataAndTime"]][ , "Time"] - 
                                  parameters[["categoricalDataAndTime"]][1, "Time"])
        return(deltaTime)
      })
        
      
      #Text display of imported file
      output$mainFileAssignmentText <- renderText({
        mainFile = uploadedFile$main
        createColoredTextForBooleanViewing(!is.null(mainFile), stateText = "Raw file name: ", invalidStateText = "Unassigned", validStateText = mainFile$name)
      })
      output$standardFileAssignmentText <- renderText({
        standardFile = uploadedFile$standard
        createColoredTextForBooleanViewing(!is.null(standardFile), stateText = "Standard file name: ", invalidStateText = "Unassigned", validStateText = standardFile$name)
      })
      output$internalStandardFileAssignmentText <- renderText({
        internalStandardFile = uploadedFile$internalStandard
        createColoredTextForBooleanViewing(!is.null(internalStandardFile), stateText = "ISTD file name: ", invalidStateText = "Unassigned", validStateText = internalStandardFile$name)
      })
      
      #Table rendering of current uploaded file
      output$uploadedFilePreviewTable <- renderTable({
        req(input$uploadedFile)
        tablePreview = read.table(input$uploadedFile$datapath, sep =input$csvDelimitation, header = FALSE)
        rowNumber <- min(nrow(tablePreview), input$fileUpload_nrow) 
        columnNumber <- min(ncol(tablePreview), input$fileUpload_ncolumn) 
        tablePreview[1:rowNumber, 1:columnNumber]
      })
      
      #Buttons to set raw, std and ISTD files
      observeEvent(input$setAsMainFile, {
        req(input$uploadedFile)
        uploadedFile$main <- input$uploadedFile
      })
      observeEvent(input$setAsStandardFile, {
        req(input$uploadedFile)
        uploadedFile$standard <- input$uploadedFile
      })
      observeEvent(input$setAsInternalStandardFile, {
        req(input$uploadedFile)
        uploadedFile$internalStandard <- input$uploadedFile
      })
      
      #Button to download the list of analyte and ISTD
      output$downloadISTDTemplate <- downloadHandler(
        filename = "ISTD_Template.csv",
        content = function(file) {
          
          req(isExtractionReady() & !is.null(!is.null(extracted$internalStandardToAnalyteAssignment)))
          
          write.table(createISTDtemplate(dataFileName = (uploadedFile$main)$datapath,
                                         sep = input$csvDelimitation),
                      file, sep = input$csvDelimitation, quote = FALSE,
                      row.names = FALSE, col.names = TRUE)
        })
      
      
      isExtractionReady <- reactive({
        !is.null(uploadedFile$main) & !is.null(uploadedFile$standard)
      })
      
      #Button to extract important information and signal of the dataframe
      observeEvent(input$extract, {
        
        req(isExtractionReady())
        
        mainFileDatapath = uploadedFile$main$datapath
        standardFileDatapath = uploadedFile$standard$datapath
        internalStandardFileDatapath = uploadedFile$internalStandard$datapath
        
        extracted$firstRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation)
        extracted$firstRowOfMain <- fillEmptyStrings(extracted$firstRowOfMain)
        
        extracted$secondRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation, skip=1)
        
        requiredStringsInFirstRowOfMain <- "Sample"
        requiredStringsInSecondRowOfMain <- c("Acq. Date-Time", "Sample Name", "Type", "Level", "CPS", "CPS RSD")
        
        if (!all(requiredStringsInFirstRowOfMain %in% extracted$firstRowOfMain) || !all(requiredStringsInSecondRowOfMain %in% extracted$secondRowOfMain)) {
          shinyalert("Impossible to extract", paste("Missing either the Sample column in the first row of the main file, or  one of the following columns in the second row: ", paste(requiredStringsInSecondRowOfMain, collapse = ', '), sep = " "), type = "error")
          return(NULL)
        }
        
        extracted$main <- read.table(mainFileDatapath, skip = 2, header = FALSE, sep = input$csvDelimitation, stringsAsFactors=FALSE)
        colnames(extracted$main) <- extracted$firstRowOfMain
        
        extracted$standard <- createStandardDataFrameFromFile(dataPath = standardFileDatapath, sep = input$csvDelimitation)
        if (is.null(extracted$standard)) return(NULL)
        
        parameters$sampleNumber <- nrow(extracted$main)
        parameters$analyteNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & !grepl("ISTD", extracted$firstRowOfMain)]
        parameters$analyteNumber <- length(parameters$analyteNames)
        parameters$internalStandardNames <- extracted$firstRowOfMain[extracted$secondRowOfMain == "CPS" & grepl("ISTD", extracted$firstRowOfMain)]
        parameters$internalStandardNumber <- length(parameters$internalStandardNames)
        
        sampleTime <- as.POSIXct(extracted$main[ , which(extracted$secondRowOfMain == "Acq. Date-Time")], format= input$dateFormat)
        if(any(is.na(sampleTime))) {
          shinyalert("Impossible to extract", paste0("Error, non-valid (NA) dates at sample number(s) ", paste(which(is.na(sampleTime)), collapse = ", "),  ". Check format in csv file and in File upload and parameters tab"), type = "error")
          return(NULL)
        }
        
        parameters[["categoricalDataAndTime"]] <- data.frame(sampleTime,
                                                             extracted$main[ , which(extracted$secondRowOfMain == "Sample Name")],
                                                             extracted$main[ , which(extracted$secondRowOfMain == "Type")],
                                                             extracted$main[ , which(extracted$secondRowOfMain == "Level")],
                                                             stringsAsFactors = FALSE
        )
        colnames(parameters[["categoricalDataAndTime"]]) <- c("Time", "Sample Name", "Type", "Level")
        
        numericalColumns <- seq(from = max(which(extracted$firstRowOfMain == "Sample")) + 1, to = ncol(extracted$main))
        extracted$main[ , numericalColumns] <- sapply(extracted$main[ , numericalColumns], as.numeric)
        
        # if (!is.null(internalStandardFileDatapath)) 
        # {
        #   extracted$internalStandardToAnalyteAssignment  <- read.table(internalStandardFileDatapath, header = TRUE, sep= input$csvDelimitation, stringsAsFactors = FALSE)
        #   process$analyteToIstdBlankRatio <- process$analyteCountsPerSecondEudc()$divideBy(process$internalStandardEudcAdaptedToAnalytes())
        # }
        # else 
        # {
        #   process$analyteToIstdBlankRatio <- process$analyteCountsPerSecondEudc()
        # }
        
        applicationState$isExtractionSuccessful <- TRUE
      })
      
      #Here we render warnings texts to help the user
      output$extractionReadyText <- renderText({
        createColoredTextForBooleanViewing(isExtractionReady(), stateText = "Data extraction ", invalidStateText = "impossible", validStateText = "ready")
      })
      
      return(
        list(extracted = extracted,
             parameters = parameters,
             applicationState = applicationState)
        )
    }
  )
}
