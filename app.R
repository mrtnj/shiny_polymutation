library(shiny)
library(ggplot2)
library(plyr)
library(reshape2)

source("simulation_functions.R")

ui <- fluidPage(
  tags$h1("Selection, mutation and drift at multiple loci"),
  column(4,
         tabsetPanel(
             tabPanel("Multiplicative",
                 sliderInput(inputId = "log_mu",
                             label = "log10(Mutation rate)",
                             value = -4, min = -6, max = -1),
                 sliderInput(inputId = "log_s",
                             label = "log10(Selection coefficient)",
                             value = -1, min = -6, max = 0),
                 sliderInput(inputId = "h",
                             label = "Dominance coefficient",
                             value = 0.25, min = 0, max = 1),
                 sliderInput(inputId = "log_loci",
                             label = "log10(Number of loci)",
                             value = 2, min = 1, max = 3),
                 actionButton(inputId = "run_button",
                             label = "Run")),
             tabPanel("Epistatic")),
         tags$p("See the code on ",
                tags$a(href = "https://github.com/mrtnj/shiny_polymutation", "GitHub"))),
  column(8,
         plotOutput(outputId = "plot_fitness"),
         plotOutput(outputId = "plot_q"))
)

server <- function(input, output) {
  simulated_data <- reactiveValues(data = NULL)

  observeEvent(input$run_button, {
    N <- 500
    s <- 10^input$log_s
    h <- input$h
    loci <- 10^input$log_loci
    mu <- 10^input$log_mu

    update_progress <- function(value) {
      incProgress(value)
    }

    withProgress(message = "Simulating", {
      sim <- sim_variation(N = N,
                           mu = mu,
                           s = s,
                           h = h,
                           loci = loci,
                           gen = 200,
                           progress_function = update_progress)
    })

    plots <- plot_simulations(sim)

    simulated_data$data <- plots
  })
  output$plot_fitness <- renderPlot({
    if (! is.null(simulated_data$data))
      print(simulated_data$data[[1]])
  })
  output$plot_q <- renderPlot({
    if (! is.null(simulated_data$data))
      print(simulated_data$data[[2]])
  })
}

# Run the application
shinyApp(ui = ui, server = server)

