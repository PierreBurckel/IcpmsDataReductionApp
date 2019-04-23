library(shiny)

source('ICPMS_ui.R')
source('ICPMS_server.R')

app <- shinyApp(ui = ICPMS_ui, server = ICPMS_server)

runApp(app, port = 4856)

