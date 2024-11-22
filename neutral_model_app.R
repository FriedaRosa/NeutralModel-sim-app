### libraries ###

library(shiny)
library(shinyWidgets)
library(spatstat)
library(terra)



### functions ###


##### NEUTRAL TEMPORAL MODEL #######
simulate_neutral_model_with_births <- function(
    n_individuals_initial,       # Number of initial individuals
    n_steps,                     # Number of timesteps
    dispersal_distance,          # Maximum dispersal distance
    grid_size,                   # Size of the grid
    cluster_sd,                  # Standard deviation for clustering
    n_clusters,                  # Number of clusters
    n_new_individuals_per_step,  # Number of new individuals added at each timestep
    dispersal_probability        # Probability of dispersal for each individual
) {
  ###### CREATE THE INITIAL POINT PATTERN ########
  # Generate cluster centers
  cluster_centers <- data.frame(
    x = runif(n_clusters, 0, grid_size),
    y = runif(n_clusters, 0, grid_size)
  )
  
  # Adjust individuals per cluster
  individuals_per_cluster <- rep(1:n_clusters, length.out = n_individuals_initial)
  
  # Generate initial locations around cluster centers
  locations <- data.frame(
    x = cluster_centers$x[individuals_per_cluster] + rnorm(n_individuals_initial, mean = 0, sd = cluster_sd),
    y = cluster_centers$y[individuals_per_cluster] + rnorm(n_individuals_initial, mean = 0, sd = cluster_sd)
  )
  
  # Ensure all points remain within the grid boundaries
  locations$x <- pmin(pmax(locations$x, 0), grid_size)
  locations$y <- pmin(pmax(locations$y, 0), grid_size)
  
  # Convert to point pattern
  pp <- ppp(locations$x, locations$y, c(0, grid_size), c(0, grid_size))
  
  # Initialize the results list
  results <- list()
  results[[1]] <- pp  # Store the initial pattern
  
  # Simulate over time
  for (step in 2:(n_steps + 1)) {
    # Determine dispersal for each individual
    dispersal_mask <- runif(length(pp$x)) < dispersal_probability
    new_x <- pp$x
    new_y <- pp$y
    
    # Apply dispersal only for individuals with dispersal_mask = TRUE
    new_x[dispersal_mask] <- new_x[dispersal_mask] + runif(sum(dispersal_mask), -dispersal_distance, dispersal_distance)
    new_y[dispersal_mask] <- new_y[dispersal_mask] + runif(sum(dispersal_mask), -dispersal_distance, dispersal_distance)
    
    # Keep points within the grid
    new_x <- pmin(pmax(new_x, 0), grid_size)
    new_y <- pmin(pmax(new_y, 0), grid_size)
    
    # Generate new individuals (births) with uniform dispersal around existing points
    if (n_new_individuals_per_step > 0) {
      birth_indices <- sample(seq_along(pp$x), n_new_individuals_per_step, replace = TRUE)
      
      # Generate random angles and distances for dispersal
      angles <- runif(n_new_individuals_per_step, 0, 2 * pi)
      distances <- runif(n_new_individuals_per_step, 0, dispersal_distance)
      
      # Calculate new positions for births
      birth_x <- pp$x[birth_indices] + distances * cos(angles)
      birth_y <- pp$y[birth_indices] + distances * sin(angles)
      
      # Keep new points within the grid
      birth_x <- pmin(pmax(birth_x, 0), grid_size)
      birth_y <- pmin(pmax(birth_y, 0), grid_size)
      
      # Combine old and new individuals
      combined_x <- c(new_x, birth_x)
      combined_y <- c(new_y, birth_y)
    } else {
      combined_x <- new_x
      combined_y <- new_y
    }
    
    # Update the point pattern
    pp <- ppp(combined_x, combined_y, c(0, grid_size), c(0, grid_size))
    
    # Store the current time step result
    results[[step]] <- pp
  }
  
  return(results)
}

############ AGGREGATIONS ON THE GRID #######################

# Function to aggregate spatial data while preserving the original extent
aggregate_spatial <- function(result_list, grid_size, scales) {
  n_steps <- length(result_list)  # Number of time steps
  aggregated_results <- list()   # Store results for each scale
  
  for (scale in scales) {
    aggregated_time_steps <- list()
    
    # Check if the scale is valid
    if (grid_size %% scale != 0) {
      stop(paste("grid_size must be divisible by scale. For grid_size =", grid_size, ", scale =", scale, "is invalid."))
    }
    
    # Loop through each time step
    for (t in seq_len(n_steps)) {
      pp <- result_list[[t]]  # Extract point pattern for the time step
      
      # Create a raster of counts at the initial resolution
      r <- rast(xmin = 0, xmax = grid_size, ymin = 0, ymax = grid_size, resolution = 1)
      
      # Rasterize the point pattern (count of points in each cell)
      counts <- rasterize(cbind(pp$x, pp$y), r, fun = max, background = 0)
      
      # Aggregate the raster using the specified scale
      counts_agg <- aggregate(counts, fact = scale, fun = max, na.rm = TRUE)
      
      aggregated_time_steps[[t]] <- counts_agg  # Store aggregated raster for the time step
    }
    
    aggregated_results[[paste0(scale, "x", scale)]] <- aggregated_time_steps
  }
  
  return(aggregated_results)
}



### shiny ##
library(shiny)

ui <- fluidPage(
  titlePanel("Neutral Temporal Model Simulation"),
  sidebarLayout(
    sidebarPanel(
      # Existing inputs
      sliderInput("n_individuals_initial", "Initial Number of Individuals:",
                  min = 1, max = 100, value = 1),
      sliderInput("n_steps", "Number of Time Steps:",
                  min = 1, max = 100, value = 50),
      sliderInput("dispersal_distance", "Dispersal Distance:",
                  min = 1, max = 100, value = 50),
      sliderInput("grid_size", "Grid Size:",
                  min = 10, max = 512, value = 256, step = 1),
      sliderInput("cluster_sd", "Cluster Standard Deviation:",
                  min = 0.1, max = 50, value = 3, step = 0.1),
      sliderInput("n_clusters", "Number of Clusters:",
                  min = 1, max = 10, value = 1),
      sliderInput("n_new_individuals_per_step", "New Individuals per Step:",
                  min = 0, max = 10, value = 1),
      sliderInput("dispersal_probability", "Dispersal Probability:",
                  min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("scaleIndex", "Aggregation Scale:",
                  min = 1, max = 9, value = 1, step = 1),
      # New Time Step slider
      sliderInput("timeStep", "Select Time Step:",
                  min = 1, max = 1, value = 1, step = 1),
      actionButton("runSim", "Run Simulation")
    ),
    mainPanel(
      plotOutput("simulationPlot"),
      plotOutput("aggregationPlot")
    )
  )
)





server <- function(input, output, session) {
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



#### RUN THE APP #####


shinyApp(ui = ui, server = server)
