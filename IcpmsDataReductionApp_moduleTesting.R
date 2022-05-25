library(tidyverse)
library(shiny)
library(shinyalert)
library(shinyjs)
library(plotly)
library(DT)
library(rstudioapi)
library(R6)
library(logging)

# install.packages("tidyverse")
# install.packages("shiny")
# install.packages("shinyalert")
# install.packages("shinyjs")
# install.packages("plotly")
# install.packages("DT")
# install.packages("rstudioapi")
# install.packages("R6")
# install.packages("logging")

moduleName <- "fileUpload"

setwd(dirname(getActiveDocumentContext()$path))  

source('IcpmsDataReductionApp_functions.R')
source('IcpmsDataReductionApp_global.R')
source(paste0('IcpmsDataReductionApp_', moduleName, '.R'))

options(shiny.error = function() { 
  logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })

ICPMS_ui <- shinyUI({
  fluidPage(
    useShinyalert(),
    shinyjs::useShinyjs(),
    do.call(paste0(moduleName, "_ui"), list())
  )
})

ICPMS_server <- function(input, output, session) {
  do.call(paste0(moduleName, "_server"), list(input = input, output = output, session = session))
}

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app,  port = 4856, launch.browser = TRUE)
