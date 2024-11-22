library(shiny)
library(shinyWidgets)
library(spatstat)
library(terra)

source("functions.R")

fluidPage(
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