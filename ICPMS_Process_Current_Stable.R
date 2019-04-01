library(shiny)
library(DT)
library(stringr)

#Functions declaration

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
  std.raw <- read.table(stdFileName, header = TRUE, sep=';', stringsAsFactors=FALSE)
  
  #Remove potential duplicates of element names in the standard dataframe
  std.raw <- removeDuplicateLines(std.raw)
  
  return(list(raw=dat.raw, std=std.raw, header_1=header_1, header_2=header_2))
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


processData <- function(p.ratio.1, p.RSD.1, e_names, std.raw, drift_ind, levelColumn, timeColumn, dtimeColumn, e_drift_choice){
  e_pt <- "[ ]{2}[A-Z]{1}[a-z]*[ ]{2}"
  e_nb <- length(e_names)
  row.names(std.raw) <- std.raw[,1]
  std.raw <- std.raw[,2:length(std.raw)]
  for (i in 1:e_nb){
    e_fn <- e_names[i]
    e_n = gsub(" ", "", str_extract(e_fn,e_pt), fixed = TRUE)
    std_c <- t(std.raw[e_n,])
    if (sum(std_c, na.rm=TRUE) != 0)
    {
      std_c <- std_c[!is.na(std_c),1, drop = FALSE]
      e_std_n.ind <- sapply(paste("^", row.names(std_c), "$", sep=""),
                            grep, make.names(levelColumn))
      e_sig <- p.ratio.1[e_std_n.ind,i] #plante
      e_drift_n.ind <- drift_ind[which(timeColumn[drift_ind] >= timeColumn[min(e_std_n.ind)])]
      e_t0 <- dtimeColumn[min(e_drift_n.ind)]
      e_dt <- dtimeColumn - e_t0
      e_dt[e_dt < 0] <- NA
      dt <- e_dt[e_drift_n.ind]
      cal_m <- lm(std_c ~ 0+e_sig)
      if (all(is.na(p.ratio.1[e_drift_n.ind,i]))){
        pred = rep(NA, nrow(p.ratio.1))
      }
      else{
        m2 <- lm(p.ratio.1[e_drift_n.ind,i] ~ poly(dt, degree=2, raw=TRUE))
        pred=predict(m2, newdata = data.frame(dt=e_dt))
      }
      if (i == 1){
        dat.drift <- data.frame(pred)
        dat.coef <- summary(cal_m)$coefficients[1,1]
      } else{
        dat.drift[i] <- pred
        dat.coef[i] <- summary(cal_m)$coefficients[1,1]
      }
      if (e_drift_choice[i]){
        dat.drift[i] <- dat.drift[i] / summary(m2)$coefficients[1,1]
        dat.drift[is.na(dat.drift[,i]), i] <- 1
      } else {dat.drift[i] <- rep(1, nrow(p.ratio.1))}
    }
    else{
      if (i == 1){
        dat.drift <- data.frame(rep(NA, nrow(p.ratio.1)))
        dat.coef <- NA
      }
      else{
        dat.drift[i] <- rep(NA, nrow(p.ratio.1))
        dat.coef[i] <- NA
      }
    }
  }
  p.ratio.2 <- p.ratio.1 / dat.drift
  p.RSD.2 <- p.RSD.1
  smp_conc <- t(t(p.ratio.2)*dat.coef)
  smp_conc[smp_conc < 0] <- "<blk"
  
  smp_RSD <- p.RSD.2
  smp_RSD[smp_conc < 0] <- "N/A"
  
  return(list(smp_conc,smp_RSD))
}





ui <- fluidPage(
  navbarPage("ICP-MS processing",
             tabPanel(
               "File upload and parameters",
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file", "File upload"),
                   actionButton("setAsRaw", "Set as raw file"),
                   actionButton("setAsStd", "Set as standard file"),
                   downloadButton("downloadISTDTemplate", "Download ISTD template"),
                   actionButton("setAsISTD", "Set as ISTD file"),
                   actionButton("extract", "Extract"),
                   htmlOutput("extract_ready_txt"),
                   htmlOutput("ISTD_not_extracted_txt")),
                 mainPanel(
                   htmlOutput("raw_assignment_txt"),
                   htmlOutput("std_assignment_txt"),
                   htmlOutput("ISTD_assignment_txt"),
                   div(style="display:inline-block",numericInput(inputId="fileUpload_nrow", label="Rows", value = 10, min = 1)),
                   div(style="display:inline-block",numericInput(inputId="fileUpload_ncolumn", label="Columns", value = 6, min = 1)),
                   tableOutput("table")
                 )
               )
             ),
             tabPanel("Index creation",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("searchIndexwhere", "In:", c("Sample Names" = "smp", "Levels" = "lvl", "Type" = "type")),
                          textInput("searchIndexwhat", "Search:", ""),
                          selectInput("searchIndexhow", "Using:", c("Regular expression" = "regexp", "Exact match" = "ematch")),
                          textInput("searchIndexName", "Index Name:", ""),
                          actionButton("indexSelectAll", "Select All"),
                          actionButton("searchIndexCreate", "Create new index")),
                        mainPanel(
                          selectInput("searchIndexDisplay", "Show:", c("Internal standards" = "ISTD", "Analytes" = "analytes")),
                          DT::DTOutput("indexTable")
                        )
                      )
             ),
             tabPanel("ISTD verification/processing",
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::useShinyjs(),
                          sliderInput("ISTDrowSlider", label = h4("Row/Column range"), min = 1, 
                                      max = 100, value = c(1, 1)),
                          sliderInput("ISTDcolSlider", label = NULL, min = 1, 
                                    max = 100, value = c(1, 1)),
                          selectInput("ISTDinteractionMode", "Select a mode:",
                                      c("View" = "view", "Process" = "process")),
                          selectInput("indexISTDchoiceWhat", "Replace what:",
                                      "", selected = ""),
                          selectInput("ISTDcorrectionMethod", "Replace calibration ISTD:",
                                      c("None" = "none", "Average blanks" = "mean",
                                        "Previous blank" = "prev")),
                          selectInput("indexISTDchoiceIn", "Replace in:",
                                      "All", selected = "All"),
                          actionButton("setISTDcorrectionMethod", "Set ISTD correction")),
                        mainPanel(
                          DT::DTOutput("ISTDtable")
                        )
                      )
             ),
             tabPanel("Blank verification/processing",
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::useShinyjs(),
                          selectInput("blkInteractionMode", "Select a mode:",
                                      c("View" = "view", "Process" = "process")),
                          selectInput("indexBlkchoiceWhat", "Replace what:",
                                      "", selected = ""),
                          selectInput("blkInterpolationMethod", "Interpolate blanks:",
                                      c("None" = "none", "Average blanks" = "mean",
                                        "Previous blank" = "prev")),
                          selectInput("indexBlkchoiceIn", "Replace in:",
                                      "All", selected = "All"),
                          actionButton("setBlkInterpolationMethod", "Set blank interpolation")),
                        mainPanel(
                          div(style="display:inline-block",sliderInput("blkRowSlider", label = h3("Row range"), min = 1, 
                                                                       max = 100, value = c(1, 1))),
                          div(style="display:inline-block",sliderInput("blkColSlider", label = h3("Column range"), min = 1, 
                                                                       max = 100, value = c(1, 1))),
                          DT::DTOutput("blkTable")
                        )
                      )
             ),
             tabPanel("Drift verification/processing",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("selectDriftIndex", "Define drift index:", "", selected = ""),
                          actionButton("setAsDriftIndex", "Set drift index"),
                          selectInput("e_drift", "Choose element:", "", selected = ""),
                          numericInput("e_ind_drift", "Element number", 1),
                          actionButton("setDriftCorrection", "Set ISTD correction"),
                          textOutput("warningDrifr")),
                        mainPanel(
                          tableOutput("smpBlkCorTable"),
                          plotOutput("driftPlot"),
                          tableOutput("test")
                        )
                      )
             ),
             tabPanel("Process",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("viewConcentrationSwitch", "Display:",
                                      c("Concentrations (ppb)" = 1, "RSD (%)" = 2, "Columnwise concatenation" = 3)),
                          selectInput("viewConcentrationIndex", "View index:",
                                      "All", selected = "All"),
                          actionButton("process", "Process data"),
                          downloadButton("downloadData", "Download table")),
                        mainPanel(
                          tableOutput("conc")
                        )
                      )
             )
  )
)

server <- function(input, output, session) {
  
  #Setting up lists of reactiveValues variables
  #uploadedFile contains the information about user-uploaded files
  uploadedFile <- reactiveValues()
  #extracted$data contains the list of important data after extraction
  extracted <- reactiveValues()
  #index contains all indexes serving as masks to display/process specific lines or columns of the raw data
  index <- reactiveValues()
  process <- reactiveValues()
  #tempDataFile stores the current file loaded by the user
  tempDataFile <- reactive({input$file})
  #extractionReady is TRUE or FALSE, depending on whether extraction can be done or not (required files assigned, variables names inputed)
  extractionReady <- reactive({!is.null(uploadedFile$raw) & !is.null(uploadedFile$std)})
  isValidISTD <- reactive({(!is.null(index$IS_CPS)) & (!is.null(extracted$data[["ISTD"]]))})
  
  #Line indexes
  #index$std <- reactive({
  #index$blk <- reactive({
  #index$drift <- reactive({
  CPS <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$CPS]})
  
  RSD <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$CPS_RSD]})
  
  IS_CPS <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$IS_CPS]})
  
  IS_CPS_RSD <- reactive({
    rawData = extracted$data[["raw"]]
    rawData[,index$IS_CPS_RSD]})
  
  elementNames <- reactive({
    header_1 = extracted$data[["header_1"]]
    header_1[index$CPS]})
  
  elementNumber <- reactive({
    length(which(index$CPS))})
  
  ISNumber <- reactive({
    length(which(index$IS_CPS))})
  
  sampleNumber <- reactive({
    nrow(extracted$data[["raw"]])})
  
  #Important columns
  timeColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    as.POSIXct(rawData[,which(header_2=="Acq. Date-Time")], format="%d/%m/%Y %H:%M")})
  
  nameColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Sample Name")]})
  
  typeColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Type")]})
  
  levelColumn <- reactive({
    rawData = extracted$data[["raw"]]
    header_2 = extracted$data[["header_2"]]
    rawData[,which(header_2=="Level")]})
  
  dtimeColumn <- reactive({
    timeColumn() - timeColumn()[1]})

  liveReplaceISTDtable <- reactive({
    req(isValidISTD())
    replaceIndexWhat = index$custom[[input$indexISTDchoiceWhat]]
    replaceIndexIn = index$custom[[input$indexISTDchoiceIn]]
    replaceValues(IS_CPS(), IS_CPS_RSD(), 1:ISNumber(), input$ISTDcorrectionMethod, replaceIndexWhat, replaceIndexIn)})
  
  ISTDmatrix <- reactive({createISTDMatrix(extracted$data[["ISTD"]], process$ISTDsignal)})
  
  process$ratio <- reactive({list(signal=CPS()/ISTDmatrix(), RSD=RSD())})
  
  liveReplaceBlkTable <- reactive({
    replaceIndexWhat = index$custom[[input$indexBlkchoiceWhat]]
    replaceIndexIn = index$custom[[input$indexBlkchoiceIn]]
    replaceValues(CPS()/ISTDmatrix(), RSD(), 1:elementNumber(), input$blkInterpolationMethod, replaceIndexWhat, replaceIndexIn)})
  
  process$ratio_cor_b <- reactive({propagateUncertainty(a=process$ratio(), b=process$blk_ratio, operation="substraction")})

  yplus <- reactive({(process$ratio_cor_b()[[1]] + process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]})
  yminus <- reactive({(process$ratio_cor_b()[[1]] - process$ratio_cor_b()[[2]]/100*process$ratio_cor_b()[[1]])[index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]})
  ###########################File import and viewing###########################
  
  #Text display of imported file
  output$raw_assignment_txt <- renderText({
    raw_file = uploadedFile$raw
    renderState(!is.null(raw_file), stateTxt = "Raw file name: ", invalidStateTxt = "Unassigned", validStateTxt = raw_file$name)
  })
  output$std_assignment_txt <- renderText({
    std_file = uploadedFile$std
    renderState(!is.null(std_file), stateTxt = "Standard file name: ", invalidStateTxt = "Unassigned", validStateTxt = std_file$name)
  })
  output$ISTD_assignment_txt <- renderText({
    ISTD_file = uploadedFile$ISTD
    renderState(!is.null(ISTD_file), stateTxt = "ISTD file name: ", invalidStateTxt = "Unassigned", validStateTxt = ISTD_file$name)
  })
  
  #Table rendering of current uploaded file
  output$table <- renderTable({
    req(tempDataFile())
    tempDataTable = read.table(tempDataFile()$datapath, sep =";", header = FALSE)
    rowNumber <- min(nrow(tempDataTable), input$fileUpload_nrow) 
    columnNumber <- min(length(tempDataTable), input$fileUpload_ncolumn) 
    tempDataTable[1:rowNumber, 1:columnNumber]
  })
  
  #Buttons to set raw, std and ISTD files
  observeEvent(input$setAsRaw, {
    req(tempDataFile())
    uploadedFile$raw <- tempDataFile()
  })
  observeEvent(input$setAsStd, {
    req(tempDataFile())
    uploadedFile$std <- tempDataFile()
  })
  observeEvent(input$setAsISTD, {
    req(tempDataFile())
    uploadedFile$ISTD <- tempDataFile()
  })
  
  #Button to download the list of analyte and ISTD
  output$downloadISTDTemplate <- downloadHandler(
    filename = "ISTD_Template.csv",
    content = function(file) {
      write.csv(createISTDtemplate((uploadedFile$raw)$datapath),
                file, sep=";", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
      })

  #Button to extract important information and signal of the dataframe
  observeEvent(input$extract, {
    #Defines name space
    raw_file = uploadedFile$raw
    std_file = uploadedFile$std
    ISTD_file = uploadedFile$ISTD
    #Allows extraction if extraction condition are met (see extractionReady reactive value)
    if (!extractionReady()){
      return()
    } else {}
    
    extracted$data <-  extractData(raw_file$datapath, std_file$datapath)
    
    index$CPS <- extracted$data[["header_2"]] == "CPS" & !grepl("ISTD", extracted$data[["header_1"]])
    index$CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & !grepl("ISTD", extracted$data[["header_1"]])
    index$IS_CPS <- extracted$data[["header_2"]] == "CPS" & grepl("ISTD", extracted$data[["header_1"]])
    index$IS_CPS_RSD <- extracted$data[["header_2"]] == "CPS RSD" & grepl("ISTD", extracted$data[["header_1"]])
    index$numericalColumns <- min(which(index$CPS)):length(extracted$data[["raw"]])
    
    extracted$data[["raw"]][,index$numericalColumns] <- sapply(extracted$data[["raw"]][,index$numericalColumns], as.numeric)
    
    ##If there exists an ISTD file that has been uploaded, extract and store it
    if (!is.null(ISTD_file)) {
      extracted$data[["ISTD"]]  <- read.table(ISTD_file$datapath, header = TRUE, sep=';', stringsAsFactors=FALSE)
    } else {}
    
    ##Defines an a priori ISTD variable that will be possible to modified in the next pane
    if (isValidISTD()) {
      process$ISTDsignal <- IS_CPS()
      process$blk_ratio <- list(signal=CPS()/ISTDmatrix(),RSD=RSD())
    } else {
      process$blk_ratio <- list(signal=CPS(),RSD=RSD())
    }
    
    index$custom[["All"]] <- which(rep(x= TRUE, sampleNumber()))
    
    ##Creates a boolean vector containing the decision of whether or not to correct signal drift for each element
    process$driftCorrectedElements <- rep(FALSE,elementNumber())
  })
  
  #Here we render warnings texts to help the user
  output$extract_ready_txt <- renderText({
    renderState(extractionReady(), stateTxt = "Data extraction ", invalidStateTxt = "impossible", validStateTxt = "ready")
  })
  output$ISTD_not_extracted_txt <- renderText({
    renderState(!(!is.null(uploadedFile$ISTD) & !is.null(extracted$data) & (isValidISTD() == FALSE)), stateTxt = "", invalidStateTxt = "Caution, ISTD not extracted", validStateTxt = NULL)
  })
  
  ##################################Index creation######################
  output$indexTable <- DT::renderDT({

    if (is.null(CPS())){return()}
    
    searchWhere = input$searchIndexwhere
    firstColumnName = ""
    
    if (searchWhere == 'smp'){
      headerRows = nameColumn()
      firstColumnName = "Sample Names"}
    else if (searchWhere == 'lvl'){
      headerRows = levelColumn()
      firstColumnName = "Levels"}
    else if (searchWhere == 'type'){
      headerRows = typeColumn()
      firstColumnName = "Type"}
    
    searchWhat = input$searchIndexwhat
    searchType = input$searchIndexhow
    displayWhat = input$searchIndexDisplay
    
    if (searchType == 'ematch'){searchWhat = paste("^", searchWhat, "$", sep="")}
    
    if (searchWhat == ""){
      index$temp <- c(1:sampleNumber())
    }
    else{
      index$temp <- grep(searchWhat, headerRows)
    }
    
    if (displayWhat == "ISTD" & !is.null(IS_CPS())){
      indexTable = cbind(headerRows[index$temp], IS_CPS()[index$temp,])
    }
    else if (displayWhat == "analytes" & !is.null(CPS())){
      indexTable = cbind(headerRows[index$temp], CPS()[index$temp,])
    }
    else {return()}
    
    names(indexTable) <- c(firstColumnName, names(indexTable)[2:length(indexTable)])
    
    indexTable
    
  }, options = list(dom = '', pageLength = sampleNumber(), ordering=T))
  
  observeEvent(input$searchIndexwhat, {
    updateTextInput(session, "searchIndexName", value = input$searchIndexwhat)
  })
  
  observeEvent(input$searchIndexCreate, {
    indexName = input$searchIndexName
    customNumIndex = c(1:sampleNumber())[index$temp][input$indexTable_rows_selected]
    if (is.null(index$custom)) {
      index$custom <- list()
    }
    index$custom[[input$searchIndexName]] <- customNumIndex
  })
  
  indexTableProxy <- DT::dataTableProxy("indexTable", session = session)
  
  observeEvent(input$indexSelectAll, {
    if (is.null(input$indexTable_rows_selected)){
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
    else if (input$indexTable_rows_selected == c(1:length(index$temp))){
      DT::selectRows(indexTableProxy, NULL)
    }
    else{
      DT::selectRows(indexTableProxy, c(1:length(index$temp)))
    }
  })
  ##################################ISTD verif/process######################
  
  observeEvent(index$custom, {
    if (is.null(index$custom)){return()}
    updateSelectInput(session,"indexISTDchoiceWhat",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceWhat",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexISTDchoiceIn",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"indexBlkchoiceIn",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"selectDriftIndex",choices=names(index$custom),names(index$custom)[1])
    updateSelectInput(session,"viewConcentrationIndex",choices=c("All", names(index$custom)),"All")
  })
  
  #Render ISTD table if all conditions are met
  output$ISTDtable <- DT::renderDT({
    if (!isValidISTD()){
      print("No available ISTD")
      return()
      }
    ISTDmode = input$ISTDinteractionMode
    lr = input$ISTDrowSlider[1]
    ur = input$ISTDrowSlider[2]
    lc = input$ISTDcolSlider[1]
    uc = input$ISTDcolSlider[2]
    
    df_view <- cbind(nameColumn(), process$ISTDsignal)
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    df_process <- cbind(nameColumn()[lr:ur], liveReplaceISTDtable()[[1]][lr:ur, lc:uc, drop = FALSE])
    names(df_process) <- c("Sample Name", names(df_process)[2:length(df_process)])
    df_process <- format(df_process, digits = 3, scientific=T)
    
    if (ISTDmode == "view") {df_view}
    else if (ISTDmode == "process") {df_process}
  }, options = list(dom = '', pageLength = sampleNumber(), ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '150px', targets = "_all"))))
  
  #Defines a proxy for changing selections in the table
  ISTDtableProxy <- DT::dataTableProxy("ISTDtable", session = session)
  
  observeEvent(input$ISTDrowSlider, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$ISTDcolSlider, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexISTDchoiceWhat, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexISTDchoiceIn, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$ISTDinteractionMode, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$ISTDcorrectionMethod, {
    if (!is.null(input$indexISTDchoiceWhat) & isValidISTD()){
      DT::selectRows(ISTDtableProxy, index$custom[[input$indexISTDchoiceWhat]])
    }
    else {}
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setISTDcorrectionMethod, {
    repIndex = index$custom[[input$indexISTDchoiceWhat]]
    if (is.null(liveReplaceISTDtable()) | is.null(repIndex) | !isValidISTD()) {return()}
    lr = min(which(input$ISTDrowSlider[1] <= repIndex))
    ur = max(which(input$ISTDrowSlider[2] >= repIndex))
    lc = input$ISTDcolSlider[1]
    uc = input$ISTDcolSlider[2]
    process$ISTDsignal[repIndex,][lr:ur, lc:uc] <- liveReplaceISTDtable()[[1]][repIndex,][lr:ur, lc:uc, drop = FALSE]
  })
  
  #Updates the slider for row and col selections when modifications are made in IS_CPS(), i.e. when the extract button is hit
  observeEvent(IS_CPS(), {
    if (!isValidISTD()){return()}
    updateSliderInput(session,"ISTDrowSlider", max=nrow(IS_CPS()), value = c(1,nrow(IS_CPS())))
    updateSliderInput(session,"ISTDcolSlider", max=length(IS_CPS()), value = c(1,length(IS_CPS())))
  })
  
  observe({
    if(input$ISTDinteractionMode == "process") {
      shinyjs::enable("setISTDcorrectionMethod")
      shinyjs::enable("ISTDrowSlider")
      shinyjs::enable("ISTDcolSlider")
      shinyjs::enable("ISTDcorrectionMethod")
      shinyjs::enable("indexISTDchoiceIn")
    }
    else if(input$ISTDinteractionMode == "view") {
      shinyjs::disable("setISTDcorrectionMethod")
      shinyjs::disable("ISTDrowSlider")
      shinyjs::disable("ISTDcolSlider")
      shinyjs::disable("ISTDcorrectionMethod")
      shinyjs::disable("indexISTDchoiceIn")
    }
  })

##################################Blank verif/process######################
  
  #Render ISTD table if all conditions are met
  output$blkTable <- DT::renderDT({
    if (is.null(CPS())){return()}
    blkMode = input$blkInteractionMode
    
    lr = input$blkRowSlider[1]
    ur = input$blkRowSlider[2]
    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]
    
    df_view <- cbind(nameColumn(), process$blk_ratio[[1]])
    names(df_view) <- c("Sample Name", names(df_view)[2:length(df_view)])
    df_view <- format(df_view, digits = 3, scientific=T)
    
    df_process <- cbind(nameColumn()[lr:ur], liveReplaceBlkTable()[[1]][lr:ur, lc:uc, drop = FALSE])
    names(df_process) <- c("Sample Name", names(df_process)[2:length(df_process)])
    df_process <- format(df_process, digits = 3, scientific=T)

    if (blkMode == "view") {df_view}
    else if (blkMode == "process") {df_process}
  }, options = list(dom = '', pageLength = sampleNumber(), ordering=F, autoWidth = TRUE, scrollX=T, columnDefs = list(list(width = '120px', targets = "_all"))))
  
  #Defines a proxy for changing selections in the table
  blkTableProxy <- DT::dataTableProxy("blkTable", session = session)
  
  observeEvent(input$blkRowSlider, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkColSlider, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceWhat, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$indexBlkchoiceIn, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInteractionMode, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  observeEvent(input$blkInterpolationMethod, {
    if (!is.null(input$indexBlkchoiceWhat)){
      DT::selectRows(blkTableProxy, index$custom[[input$indexBlkchoiceWhat]])
    }
    else {}
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$setBlkInterpolationMethod, {
    repIndex = index$custom[[input$indexBlkchoiceWhat]]
    if (is.null(liveReplaceBlkTable())| is.null(repIndex)) {return()}
    lr = min(which(input$blkRowSlider[1] <= repIndex))
    ur = max(which(input$blkRowSlider[2] >= repIndex))
    lc = input$blkColSlider[1]
    uc = input$blkColSlider[2]
    process$blk_ratio[[1]][repIndex,][lr:ur, lc:uc] <- liveReplaceBlkTable()[[1]][repIndex,][lr:ur, lc:uc, drop = FALSE]
    process$blk_ratio[[2]][repIndex,][lr:ur, lc:uc] <- liveReplaceBlkTable()[[2]][repIndex,][lr:ur, lc:uc, drop = FALSE]
  })
  
  #Updates the slider for row and col selections when modifications are made in IS_CPS(), i.e. when the extract button is hit
  observeEvent(CPS(), {
    if (is.null(CPS())){return()}
    updateSliderInput(session,"blkRowSlider", max=sampleNumber(), value = c(1,sampleNumber()))
    updateSliderInput(session,"blkColSlider", max=elementNumber(), value = c(1,elementNumber()))
  })
  
  observe({
    if(input$blkInteractionMode == "process") {
      shinyjs::enable("setBlkInterpolationMethod")
      shinyjs::enable("blkRowSlider")
      shinyjs::enable("blkColSlider")
      shinyjs::enable("blkInterpolationMethod")
      shinyjs::enable("indexBlkchoiceIn")
    }
    else if(input$blkInteractionMode == "view") {
      shinyjs::disable("setBlkInterpolationMethod")
      shinyjs::disable("blkRowSlider")
      shinyjs::disable("blkColSlider")
      shinyjs::disable("blkInterpolationMethod")
      shinyjs::disable("indexBlkchoiceIn")
    }
  })
  
  #Assigns the current state of the ISTD table to the ISTD variable that will be used for calculations
  observeEvent(input$searchBlkSelect, {
    DT::selectRows(liveReplaceBlkTableProxy, NULL)
  })
  
  ##################################Drift verif/process######################
  observeEvent(input$setAsDriftIndex, {
    if (input$selectDriftIndex != ""){
      index$drift <- index$custom[[input$selectDriftIndex]]
    }
  })
  
  output$smpBlkCorTable <- renderTable({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(index$drift) | !is.integer(input$e_ind_drift)){return()}
    driftSignal <- process$ratio_cor_b()[[1]]
    if (input$e_ind_drift < 5){
      lc <- 1
    }
    else{
      lc <- max(which((1:input$e_ind_drift)%%5 == 0))
    }
    uc <- min(lc + 4, elementNumber())
    driftTable <- cbind(nameColumn()[index$drift],driftSignal[index$drift,lc:uc])
    names(driftTable) <- c("Sample Name", names(driftTable)[2:length(driftTable)])
    driftTable
    
  }, digits = -2)
  
  observeEvent(elementNames(), {
    if (is.null(elementNames())){return()}
    updateSelectInput(session,"e_drift",choices=elementNames(),selected=elementNames()[1])
  })
  
  observeEvent(input$e_drift, {
    if (is.null(input$e_drift)){return()}
    updateNumericInput(session, "e_ind_drift", value = grep(input$e_drift,elementNames(),fixed=TRUE))
  })
  
  observeEvent(input$e_ind_drift, {
    if (is.null(input$e_ind_drift)){return()}
    updateSelectInput(session,"e_drift",selected=elementNames()[input$e_ind_drift])
  })
  
  output$driftPlot <- renderPlot({
    if (is.null(process$ratio_cor_b()[[1]]) | is.null(index$drift)){return()}
    driftTime = dtimeColumn()[index$drift]
    driftValue = process$ratio_cor_b()[[1]][index$drift,grep(input$e_drift, elementNames(),fixed=TRUE)]
    plot(x=driftTime,y=driftValue,xlab = "Time (seconds)", ylab = input$e_drift, pch=21, bg="lightblue")

    if (process$driftCorrectedElements[input$e_ind_drift] == TRUE){
      time_0 = driftTime[1]
      time_f = driftTime[length(driftTime)]
      timeInterval = seq(from=as.numeric(time_0), to=as.numeric(time_f), by=as.numeric((time_f - time_0)/100))
      
      dt = as.numeric(driftTime)
      m2 <- lm(driftValue ~ poly(dt, degree=2, raw=TRUE))
      
      pred=predict(m2, newdata = data.frame(dt = timeInterval))
      
      lines(x=timeInterval, y=pred, col="red")
    }
    arrows(dtimeColumn()[index$drift], yminus(), dtimeColumn()[index$drift], yplus(), length=0.05, angle=90, code=3)
  })
  
  observeEvent(input$setDriftCorrection, {
    if (is.null(process$driftCorrectedElements)){return()}
    process$driftCorrectedElements[input$e_ind_drift] <- !process$driftCorrectedElements[input$e_ind_drift]
  })
  
  output$test <- renderTable({
    if (is.null(process$driftCorrectedElements)){return()}
    process$driftCorrectedElements
  })
  ####################Process
  observeEvent(input$process, {
    if (is.null(process$driftCorrectedElements)){return()}
    process$conc <- processData(process$ratio_cor_b()[[1]], process$ratio_cor_b()[[2]], elementNames(), extracted$data[["std"]], index$drift,levelColumn(), timeColumn(), dtimeColumn(), process$driftCorrectedElements)
  })
  
  output$conc <- renderTable({
    if (is.null(process$conc)){return()}
    
    displayWhat = strtoi(input$viewConcentrationSwitch)
    displayIndex = input$viewConcentrationIndex
    
    if (displayWhat <= 2){displayMatrix <- process$conc[[displayWhat]]}
    else if (displayWhat == 3){
      displayMatrix <- mergeMatrixes(matrix1 = process$conc[[1]], matrix2 = process$conc[[2]], name1=NULL, name2="RSD (%)")
    }
    else {}
    
    cbind(nameColumn(), displayMatrix)[index$custom[[displayIndex]],]
  })
  
  #Button to download the list of analyte and ISTD
  output$downloadData <- downloadHandler(
    filename = paste("data_", input$viewConcentrationIndex, ".csv", sep = ""),
    content = function(file) {
      selectedIndex = index$custom[[input$viewConcentrationIndex]]
      
      combinedConcRSD <- mergeMatrixes(matrix1 = process$conc[[1]], matrix2 = process$conc[[2]], name1=NULL, name2="RSD (%)")

      if (strtoi(input$viewConcentrationSwitch) <= 2){
        write.csv(cbind(nameColumn(),process$conc[[strtoi(input$viewConcentrationSwitch)]])[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
      else if (strtoi(input$viewConcentrationSwitch) == 3){
        write.csv(cbind(nameColumn(), combinedConcRSD)[selectedIndex,],
                  file, sep=";", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }
    })
  
}

app <- shinyApp(ui = ui, server = server)

runApp(app, port = 4856)

