library(tidyverse)
library(shiny)
library(shinyalert)
library(shinyjs)
library(plotly)
library(DT)
library(rstudioapi)
library(R6)

# install.packages("tidyverse")
# install.packages("shiny")
# install.packages("shinyalert")
# install.packages("shinyjs")
# install.packages("plotly")
# install.packages("DT")
# install.packages("rstudioapi")
# install.packages("R6")

setwd(dirname(getActiveDocumentContext()$path))  

source('helpers/functions.R')

source('modules/fileUpload.R')
source('modules/indexCreation.R')
source('modules/interferenceCorrection.R')
source('modules/blankProcessing.R')
source('modules/driftProcessing.R')
source('modules/results.R')

ICPMS_ui <- shinyUI({
  fluidPage(
    useShinyalert(),
    shinyjs::useShinyjs(),
    navbarPage("ICP-MS processing",
               fileUpload_ui("fileUpload"),
               indexCreation_ui("indexCreation"),
               interferenceCorrection_ui(),
               blankProcessing_ui("blankProcessing"),
               driftProcessing_ui(),
               results_ui()
    )
  )
})

ICPMS_server <- function(input, output, session) {
  fileUpload <- fileUpload_server("fileUpload")
  indexCreation <- indexCreation_server("indexCreation", fileUpload)
  blankProcessing_server("blankProcessing", fileUpload, indexCreation)
  
}

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app,  port = 4856, launch.browser = TRUE)
