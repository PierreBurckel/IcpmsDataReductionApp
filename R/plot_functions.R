# library(plotly)
# 
# categoricalDataAndTime
# data_eudc <- reactive({data_to_plot[ , input$data_to_plot]})
# data_df <- cbind(categoricalDataAndTime, data_eudc$getEstimation()) %>%
#   as.data.frame()
# plotData(data_df)
# 
# plotData <- function(data_df) {
#   element_names <- data_eudc$getElementFullNames()
#   data_df <- as.data.frame(data_eudc$getEstimation())
#   fig <- plot_ly(data, x = ~x, y = ~trace_0, name = 'trace 0', type = 'scatter', mode = 'lines+markers')
# }
# 
# trace_0 <- rnorm(100, mean = 5)
# 
# trace_1 <- rnorm(100, mean = 0)
# 
# trace_2 <- rnorm(100, mean = -5)
# 
# x <- c(1:100)
# 
# 
# data <- data.frame(x, trace_0, trace_1, trace_2)
# 
# 
# fig <- plot_ly(data, x = ~x, y = ~trace_0, name = 'trace 0', type = 'scatter', mode = 'lines') 
# 
# fig <- fig %>% add_trace(y = ~trace_1, name = 'trace 1', mode = 'lines+markers') 
# 
# fig <- fig %>% add_trace(y = ~trace_2, name = 'trace 2', mode = 'markers')
# 
# 
# fig