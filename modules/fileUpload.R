#' Module responsible for uploading data file(s) and extracting information (data 
#' tables, parameters)
#' 
#' @return A list of reactiveValues @extracted @parameters
#' 
#' extracted$ @main (data.frame) @firstRowOfMain (chr)
#' @secondRowOfMain (chr) @standard (data.frame)
#' 
#' parameters$ @deltaTime() (num) @sampleNumber (int) @analyteNames (chr) 
#' @analyteNumber (int) @internalStandardNames (chr) 
#' @internalStandardNumber (int) @categoricalDataAndTime (data.frame)

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
        htmlOutput(ns("extractedText")),
        selectInput(inputId = ns("csvDelimitation"), label = "Select csv field separator (the separator will be used during extraction,
                                 you therefore need to use the same separator for all your files)",
                    choices = c(";", ","), selected = ";"),
        textInput(inputId = ns("dateFormat"), label = "Date format", value = "%d/%m/%Y %H:%M")),
      mainPanel(
        htmlOutput(ns("mainFileAssignmentText")),
        htmlOutput(ns("standardFileAssignmentText")),
        htmlOutput(ns("internalStandardFileAssignmentText")),
        div(style="display:inline-block", numericInput(inputId=ns("fileUpload_nrow"), label="Rows", value = 10, min = 1)),
        div(style="display:inline-block", numericInput(inputId=ns("fileUpload_ncolumn"), label="Columns", value = 6, min = 1)),
        tableOutput(ns("uploadedFilePreviewTable"))
      )
    )
  )
}

fileUpload_server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      # Used within the module, stores input$file (data.frame, $name, $datapath, ...)
      # Used in a reactive context by renderText
      uploadedFile <- reactiveValues()
      
      # Returned in this module
      extracted <- reactiveValues()
      parameters <- reactiveValues()
      
      # State variables
      isExtractionReady <- reactive({
        !is.null(uploadedFile$main) & !is.null(uploadedFile$standard)
      })
      
      #Text display of imported file and 
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
      
      output$extractionReadyText <- renderText({
        createColoredTextForBooleanViewing(isExtractionReady(), stateText = "Data extraction ", invalidStateText = "impossible", validStateText = "ready")
      })
      
      output$extractedText <- renderText({
        createColoredTextForBooleanViewing(isExtracted(), stateText = "Data is ", invalidStateText = "not extracted", validStateText = "extracted")
      })
      
      #Table rendering of current uploaded file
      output$uploadedFilePreviewTable <- renderTable({ req(input$uploadedFile)
        tablePreview = read.table(input$uploadedFile$datapath, sep =input$csvDelimitation, header = FALSE)
        rowNumber <- min(nrow(tablePreview), input$fileUpload_nrow) 
        columnNumber <- min(ncol(tablePreview), input$fileUpload_ncolumn)
        tablePreview[seq(rowNumber), seq(columnNumber)]
      })
      
      #Buttons to set main, std and ISTD files
      observeEvent(input$setAsMainFile, { req(input$uploadedFile)
        uploadedFile$main <- input$uploadedFile
      })
      observeEvent(input$setAsStandardFile, { req(input$uploadedFile)
        uploadedFile$standard <- input$uploadedFile
      })
      observeEvent(input$setAsInternalStandardFile, { req(input$uploadedFile)
        uploadedFile$internalStandard <- input$uploadedFile
      })
      
      #Button to download the list of analyte and ISTD
      output$downloadISTDTemplate <- downloadHandler(
        filename = "ISTD_Template.csv",
        content = function(file) { req(isExtractionReady()) 
          write.table(createISTDtemplate(dataFileName = (uploadedFile$main)$datapath,
                                         sep = input$csvDelimitation),
                      file, sep = input$csvDelimitation, quote = FALSE,
                      row.names = FALSE, col.names = TRUE)
        })
      
      # Button to extract important information and signal of the dataframe
      # Returns false if extraction not complete and true if extraction complete
      isExtracted <- eventReactive(input$extract, { req(isExtractionReady())
        
        browser()
        
        # Datapath assignment
        mainFileDatapath = uploadedFile$main$datapath
        standardFileDatapath = uploadedFile$standard$datapath
        internalStandardFileDatapath = uploadedFile$internalStandard$datapath
        
        # Extraction of the first row of the main data.frame and filling of its empty elements
        # Contains element names and a sample tag for all categorical data
        firstRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation)
        firstRowOfMain <- fillEmptyStrings(firstRowOfMain)
        
        # Extraction of the second row of the main data.frame
        # Contains specific categorical tags and CPS and CPS RSD designations
        secondRowOfMain <- scan(mainFileDatapath, nlines = 1, what = character(),sep = input$csvDelimitation, skip=1)
        
        requiredStringsInFirstRowOfMain <- "Sample"
        requiredStringsInSecondRowOfMain <- c("Acq. Date-Time", "Sample Name", "Type", "Level", "CPS", "CPS RSD")
        
        if (!all(requiredStringsInFirstRowOfMain %in% firstRowOfMain) || !all(requiredStringsInSecondRowOfMain %in% secondRowOfMain)) {
          shinyalert("Impossible to extract", paste("Missing either the Sample column in the first row of the main file, or  one of the following columns in the second row: ", paste(requiredStringsInSecondRowOfMain, collapse = ', '), sep = " "), type = "error")
          return(FALSE)
        }
        
        # Extraction of the body of the main data.frame and column name assignment
        extracted$main <- read.table(mainFileDatapath, skip = 2, header = FALSE, sep = input$csvDelimitation, stringsAsFactors=FALSE)
        column_names <- paste(firstRowOfMain, secondRowOfMain, sep = "_")
        if (!identical(column_names, unique(column_names))) {
          shinyalert("Impossible to extract, column names are not unique. Delete the duplicated columns in the csv file and upload the new file prior to extraction.", type = "error")
          return(FALSE)
        }
        colnames(extracted$main) <- column_names
        extracted$main <- as_tibble(extracted$main)
        
        # Extraction of the standard data.frame
        extracted$standard <- createStandardDataFrameFromFile(dataPath = standardFileDatapath, sep = input$csvDelimitation)
        if (is.null(extracted$standard)) return(FALSE)
        
        # Extraction of the internal standard assignment information
        if (!is.null(internalStandardFileDatapath)) {
          extracted$internalStandardToAnalyteAssignment <- read.table(internalStandardFileDatapath,
                                                                      skip = 1, header = FALSE,
                                                                      sep = input$csvDelimitation,
                                                                      stringsAsFactors=FALSE)
        }
        
        
        
        # Parameters computation
        parameters$sampleNumber <- nrow(extracted$main)
        parameters$analyteNames <- firstRowOfMain[secondRowOfMain == "CPS" & !grepl("ISTD", firstRowOfMain)]
        parameters$analyteNumber <- length(parameters$analyteNames)
        parameters$internalStandardNames <- firstRowOfMain[secondRowOfMain == "CPS" & grepl("ISTD", firstRowOfMain)]
        parameters$internalStandardNumber <- length(parameters$internalStandardNames)
        parameters$csvDelimitation <- input$csvDelimitation
        
        sampleTime <- extracted$main %>% 
          select("Sample_Acq. Date-Time") %>% 
          mutate_all(as.POSIXct, format = input$dateFormat)
        
        if(any(is.na(sampleTime))) {
          shinyalert("Impossible to extract", paste0("Error, non-valid (NA) dates at sample number(s) ", paste(which(is.na(sampleTime)), collapse = ", "),  ". Check format in csv file and in File upload and parameters tab"), type = "error")
          return(FALSE)
        }
        
        parameters[["categoricalDataAndTime"]] <- tibble(sampleTime,
                                                         extracted$main %>% 
                                                           select(paste("Sample", c("Sample Name", "Type", "Level"), sep = "_"))
        )
        colnames(parameters[["categoricalDataAndTime"]]) <- c("Time", "Sample Name", "Type", "Level")
        
        # Conversion to numeric of CPS and CPS RSD columns
        numericalColumns <- seq(from = max(which(firstRowOfMain == "Sample")) + 1, to = ncol(extracted$main))
        extracted$main[ , numericalColumns] <- sapply(extracted$main[ , numericalColumns], as.numeric)
        
        # Assessing that the number of CPS, CPS RSD and unique element names are coherent
        CPS_number <- sum(secondRowOfMain == "CPS")
        CPS_RSD_number <- sum(secondRowOfMain == "CPS RSD")
        unique_element_number <- length(unique(firstRowOfMain[numericalColumns]))
        if (CPS_number != CPS_RSD_number | CPS_number != unique_element_number) {
          shinyalert("Impossible to extract, number of columns containing the keywords `CPS`, `CPS RSD` and unique element names in the first line are not equal. Check column names in the csv file.", type = "error")
          return(FALSE)
        }
        return(TRUE)
      })
      
      parameters$deltaTime <- reactive({
        # categoricalDataAndTime are tibbles, tibble format needs to be dropped to perform operations
        # deltaTime is a numeric vector and gives the time in seconds since first sample
        # as.numeric removes the timestamp and converts to pure numeric vector
        deltaTime <- as.numeric(parameters[["categoricalDataAndTime"]][, "Time",drop=TRUE] - 
                                  parameters[["categoricalDataAndTime"]][1, "Time",drop=TRUE])
        return(deltaTime)
      })
      
      return(
        list(extracted = extracted,
             parameters = parameters)
        )
    }
  )
}
