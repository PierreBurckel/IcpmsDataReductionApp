library(shiny)
library(shinyalert)
library(plotly)
library(DT)
library(stringr)
library(plotly)
library(shinyjs)
library(stringr)
library(rstudioapi)
library(ggplot2)
library(R6)
library(logging)

# install.packages("shiny")
# install.packages("shinyalert")
# install.packages("plotly")
# install.packages("DT")
# install.packages("stringr")
# install.packages("plotly")
# install.packages("shinyjs")
# install.packages("stringr")
# install.packages("rstudioapi")
# install.packages("ggplot2")
# install.packages("R6")

setwd(dirname(getActiveDocumentContext()$path))  

options(shiny.error = function() { 
  logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })

source('IcpmsDataReductionApp_ui.R')
source('IcpmsDataReductionApp_server.R')
source('IcpmsDataReductionApp_functions.R')

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app,  port = 4856, launch.browser = TRUE)
