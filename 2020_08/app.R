library(shiny)
library(tidyverse)

# A very nice explanation of how reactivity in Shiny works - 
# https://shiny.rstudio.com/articles/understanding-reactivity.html

df <- read_csv("df.csv")
genes <- pull(df, Gene.symbol)

ui <- fluidPage(titlePanel("Outlier Detector"),
                sidebarLayout(
                    sidebarPanel(selectInput("gene", "Gene:",
                                             choices = genes)),
                    mainPanel(textOutput("text"),
                              plotOutput("plot"))
                ))

server <- function(input, output) {
    output$plot <- renderPlot({
        df %>%
            filter(Gene.symbol == input$gene) %>%
            pivot_longer(
                cols = contains("GSM"),
                names_to = "Sample",
                values_to = "Expression"
            ) %>%
            ggplot(aes(Sample, Expression)) +
            geom_point(col = "darkblue", size = 5) +
            theme_bw()
    })
    
    output$text <- renderText({
        paste0("Expression levels for gene: ", input$gene)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
