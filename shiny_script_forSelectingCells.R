library(shiny)
library(plotly)
library(ggplot2)

p <- ggplot(metadata, aes(x, y, color = celltype)) +
  geom_point(size = 2)

ui <- fluidPage(
  plotlyOutput("plot"),
  h3("Selected cells"),
  tableOutput("selected"),
  downloadButton("download", "Download Selected Cells")
)

server <- function(input, output, session) {
  
  # Store selected cells in a reactive expression
  selected_cells <- reactive({
    d <- event_data("plotly_selected")
    if (is.null(d)) return(NULL)
    metadata[d$pointNumber + 1, ]
  })
  
  output$plot <- renderPlotly({
    ggplotly(p) %>% layout(dragmode = "lasso")
  })
  
  output$selected <- renderTable({
    selected_cells()
  })
  
  # Download handler to save selected cells to a CSV
  output$download <- downloadHandler(
    filename = function() {
      paste0("selected_cells_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(selected_cells(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)

