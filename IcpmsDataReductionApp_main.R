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

setwd(dirname(getActiveDocumentContext()$path))  

options(shiny.error = function() { 
  logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })

source('IcpmsDataReductionApp_ui.R')
source('IcpmsDataReductionApp_server.R')
source('IcpmsDataReductionApp_functions.R')

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app,  port = 4856, launch.browser = TRUE)
