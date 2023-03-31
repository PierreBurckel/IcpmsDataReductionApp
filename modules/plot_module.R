plot_module_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tableOutput(ns("module_table"))
  )
}

plot_module_server <- function(id, data_to_plot) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      ns <- session$ns
      
      output$module_table <- renderTable({
        cbind(
          names(data_to_plot()),
          t(as.data.frame(lapply(data_to_plot(),
                                 paste,
                                 collapse = ", ")
                          )
            )
        )
      })
    })
}