C_LETTER_KEYCODE <- 67

uploadedFile <- reactiveValues()
extracted <- reactiveValues()
rowIndexInMain <- reactiveValues()
process <- reactiveValues()
modifiers <- reactiveValues(blank = list(),
                            interference = list())
parameters <- reactiveValues()
applicationState <- reactiveValues(isExtractionSuccessful = FALSE)

shinyjs::disable(selector = '.navbar-nav a[data-value="Index creation"')
shinyjs::disable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
shinyjs::disable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
shinyjs::disable(selector = '.navbar-nav a[data-value="Process"')

observeEvent(applicationState$isExtractionSuccessful, {
  if (applicationState$isExtractionSuccessful == TRUE) {
    shinyjs::enable(selector = '.navbar-nav a[data-value="Index creation"')
    shinyjs::enable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
    shinyjs::enable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
    shinyjs::enable(selector = '.navbar-nav a[data-value="Process"')
  }
  else if (applicationState$isExtractionSuccessful == FALSE) {
    shinyjs::disable(selector = '.navbar-nav a[data-value="Index creation"')
    shinyjs::disable(selector = '.navbar-nav a[data-value="Blank verification/processing"')
    shinyjs::disable(selector = '.navbar-nav a[data-value="Drift verification/processing"')
    shinyjs::disable(selector = '.navbar-nav a[data-value="Process"')
  }
})