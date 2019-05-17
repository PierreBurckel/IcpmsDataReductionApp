library(shiny)

source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_ui.R')
source('C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/ICPMS_server.R')

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app, port = 4856)
