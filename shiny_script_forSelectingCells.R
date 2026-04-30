library(shiny)
library(plotly)
library(ggplot2)


ui <- fluidPage(
  # Add a source name to the plotly output
  plotlyOutput("plot"),
  hr(),
  h3("Selection Status"),
  textOutput("cell_count"),
  downloadButton("download", "Download Selected Cells (CSV)")
)

server <- function(input, output, session) {
  
  # Generate the Plotly object
  output$plot <- renderPlotly({
    p <- ggplot(metadata, aes(x = x, y = y, color = Celltype, key = Cell_ID)) +
      geom_point(size = 2) +
      theme_minimal()
    
    # 'source = "A"' tells Shiny exactly which plot to watch
    ggplotly(p, source = "A") %>% 
      layout(dragmode = "lasso")
  })
  
  # Capture the selection
  selected_df <- reactive({
    # Match the source "A" from above
    d <- event_data("plotly_selected", source = "A")
    
    # If nothing is selected, return the whole metadata or NULL
    if (is.null(d) || length(d$key) == 0) {
      return(NULL)
    }
    
    # Filter metadata based on the keys captured by the lasso
    return(metadata[metadata$Cell_ID %in% d$key, ])
  })
  
  # Show a count so you know it's working before you download
  output$cell_count <- renderText({
    if (is.null(selected_df())) {
      "No cells selected. Use the Lasso or Box tool on the plot."
    } else {
      paste(nrow(selected_df()), "cells selected.")
    }
  })
  
  # The Download Handler
  output$download <- downloadHandler(
    filename = function() {
      paste0("selected_cells_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(selected_df())
      
      final_df <- selected_df()
      rownames(final_df) <- final_df$Cell_ID
      
      write.csv(final_df, file, row.names = TRUE)
    }
  )
}

shinyApp(ui, server)
