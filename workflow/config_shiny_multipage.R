# Copy-catted from https://raw.githubusercontent.com/sdparekh/zUMIs/refs/heads/main/zUMIs-config_shiny.R
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(yaml)
library(shinyBS)
library(glue)
library(htmltools)
source("scripts/app_panels.R")
source("scripts/app_server.R")

# Define UI for application
ui <- fluidPage(
  # Application title
  theme = shinythemes::shinytheme("simplex"),
  titlePanel("PETRI-seq config: Generate configuration yaml file"),
  navlistPanel(id = "mainNav", widths = c(2, 10),
               general_option_panel,
               sample_configuration_panel,
               load_yaml_panel,
               save_yaml_panel
               )
)

# Run the application
shinyApp(ui = ui, server = server)
