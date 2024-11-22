library(shiny)
library(shinyWidgets)
library(spatstat)
library(terra)

source("functions.R")

function(input, output, session) {
  scales_vector <- c(1, 2, 4, 8, 16, 32, 64, 128, 256)
  
  # Reactive expression to run the simulation when the button is clicked
  simulationResult <- eventReactive(input$runSim, {
    # Run the simulation with user inputs
    result_list <- simulate_neutral_model_with_births(
      n_individuals_initial = input$n_individuals_initial,
      n_steps = input$n_steps,
      dispersal_distance = input$dispersal_distance,
      grid_size = input$grid_size,
      cluster_sd = input$cluster_sd,
      n_clusters = input$n_clusters,
      n_new_individuals_per_step = input$n_new_individuals_per_step,
      dispersal_probability = input$dispersal_probability
    )
    return(result_list)
  })
  
  # Reactive expression for aggregation
  aggregatedData <- reactive({
    result_list <- simulationResult()
    req(result_list)  # Ensure result_list is available
    scale_numeric <- scales_vector[input$scaleIndex]
    
    # Validate that grid_size is divisible by scale_numeric
    if (input$grid_size %% scale_numeric != 0) {
      showNotification("Grid Size must be divisible by Aggregation Scale.", type = "error")
      return(NULL)
    }
    
    scales <- c(scale_numeric)
    aggregated_data <- aggregate_spatial(result_list, input$grid_size, scales)
    return(aggregated_data)
  })
  
  # Update the timeStep slider's max and value after simulation
  observeEvent(simulationResult(), {
    n_steps <- length(simulationResult())
    updateSliderInput(session, "timeStep",
                      min = 1,
                      max = n_steps,
                      value = 1,
                      step = 1)
  })
  
  # Render the simulation plot
  output$simulationPlot <- renderPlot({
    result_list <- simulationResult()
    req(result_list)  # Ensure result_list is available
    
    time_step <- input$timeStep
    pp <- result_list[[time_step]]
    plot(pp, main = paste("Simulation at Time Step", time_step))
  })
  
  # Render the aggregation plot
  output$aggregationPlot <- renderPlot({
    aggregated_data <- aggregatedData()
    req(aggregated_data)  # Ensure data is available
    
    scale_numeric <- scales_vector[input$scaleIndex]
    scale_label <- paste0(scale_numeric, "x", scale_numeric)
    
    time_step <- input$timeStep
    aggregated_raster <- aggregated_data[[scale_label]][[time_step]]
    plot(aggregated_raster, main = paste("Aggregation Scale:", scale_label, "\nTime Step:", time_step))
  })
}

