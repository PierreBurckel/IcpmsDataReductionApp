library(GMCM)
library(tidyverse)
library(shiny)
library(shinyalert)
library(shinyjs)
library(plotly)
library(DT)
library(rstudioapi)
library(R6)

# install.packages("GMCM")
# install.packages("tidyverse")
# install.packages("shiny")
# install.packages("shinyalert")
# install.packages("shinyjs")
# install.packages("plotly")
# install.packages("DT")
# install.packages("rstudioapi")
# install.packages("R6")

# setwd(dirname(getActiveDocumentContext()$path))  

source('R/functions.R')

source('modules/fileUpload.R')
source('modules/reactiveExpressions.R')
source('modules/indexCreation.R')
source('modules/interferenceCorrection.R')
source('modules/blankProcessing.R')
source('modules/driftProcessing.R')
source('modules/results.R')
source('modules/plot_module.R')

ICPMS_ui <- shinyUI({
  fluidPage(
    useShinyalert(),
    shinyjs::useShinyjs(),
    navbarPage("ICP-MS processing",
               fileUpload_ui("fileUpload"),
               indexCreation_ui("indexCreation"),
               interferenceCorrection_ui("interferenceCorrection"),
               blankProcessing_ui("blankProcessing"),
               driftProcessing_ui("driftProcessing"),
               results_ui("results"),
               reactiveExpressions_ui("reactiveExpressions")
    )
  )
})

ICPMS_server <- function(input, output, session) {
  fileUpload <- fileUpload_server("fileUpload")
  reactiveExpressions <- reactiveExpressions_server("reactiveExpressions", fileUpload, indexCreation, interferenceCorrection, blankProcessing, driftProcessing)
  indexCreation <- indexCreation_server("indexCreation", fileUpload, reactiveExpressions)
  interferenceCorrection <- interferenceCorrection_server("interferenceCorrection", fileUpload, reactiveExpressions, indexCreation)
  blankProcessing <- blankProcessing_server("blankProcessing", fileUpload, reactiveExpressions, indexCreation)
  driftProcessing <- driftProcessing_server("driftProcessing", fileUpload, reactiveExpressions, indexCreation)
  results_server("results", indexCreation, reactiveExpressions, fileUpload)

}

shinyApp(ui = ICPMS_ui, server = ICPMS_server)

# app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)
# 
# runApp(app,  port = 4856, launch.browser = TRUE)
