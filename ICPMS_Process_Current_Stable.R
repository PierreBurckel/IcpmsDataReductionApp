library(shiny)
library(shinyalert)
library(plotly)
library(DT)
library(stringr)
library(plotly)
library(shinyjs)
library(stringr)
library(rstudioapi) 

# source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_ui.R')
# source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_server.R')
# source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_functions.R')

setwd(dirname(getActiveDocumentContext()$path))  

source('ICPMS_ui.R')
source('ICPMS_server.R')
source('ICPMS_functions.R')

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app,  port = 4856, launch.browser = TRUE)
