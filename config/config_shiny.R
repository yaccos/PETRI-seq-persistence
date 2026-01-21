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
library(purrr)

resolve_config_dir <- function() {
  is_config_dir <- function(path) {
    if (!nzchar(path)) {
      return(FALSE)
    }
    dir.exists(path) && file.exists(file.path(path, "config_shiny.R"))
  }

  ascend_to_config <- function(path) {
    if (!nzchar(path)) {
      return(NULL)
    }
    if (!file.exists(path) && !dir.exists(path)) {
      return(NULL)
    }
    current <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!dir.exists(current)) {
      current <- dirname(current)
    }
    repeat {
      if (is_config_dir(current)) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
    NULL
  }

  candidates <- character()

  trailing_args <- commandArgs(trailingOnly = TRUE)
  if (length(trailing_args) > 0) {
    arg_paths <- trailing_args[file.exists(trailing_args)]
    candidates <- c(candidates, arg_paths)
  }

  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  for (arg in full_args) {
    if (startsWith(arg, file_arg)) {
      candidates <- c(candidates, sub(file_arg, "", arg))
    }
  }

  option_candidates <- c(getOption("shiny.appdir"), getOption("shiny.fileName"))
  option_candidates <- option_candidates[!vapply(option_candidates, is.null, logical(1))]
  candidates <- c(candidates, option_candidates)

  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    frame <- frames[[idx]]
    candidate <- NULL
    if (exists("ofile", envir = frame, inherits = FALSE)) {
      candidate <- get("ofile", envir = frame, inherits = FALSE)
    } else if (exists("fileName", envir = frame, inherits = FALSE)) {
      candidate <- get("fileName", envir = frame, inherits = FALSE)
    } else {
      srcfile <- attr(frame, "srcfile")
      if (!is.null(srcfile) && !is.null(srcfile$filename)) {
        candidate <- srcfile$filename
      }
    }
    if (!is.null(candidate) && nzchar(candidate)) {
      candidates <- c(candidates, candidate)
    }
  }

  env_candidates <- c(Sys.getenv("PETRISEQ_ROOT"), Sys.getenv("PROJECT_ROOT"), Sys.getenv("SHINY_APP_PATH"))
  candidates <- c(candidates, env_candidates)

  candidates <- unique(candidates[nzchar(candidates)])

  for (candidate in candidates) {
    config_dir <- ascend_to_config(candidate)
    if (!is.null(config_dir)) {
      return(config_dir)
    }
  }

  fallback <- ascend_to_config(getwd())
  if (!is.null(fallback)) {
    return(fallback)
  }

  stop("Unable to locate config directory containing config_shiny.R. Checked candidates: ",
       paste(candidates, collapse = ", "))
}

config_dir <- resolve_config_dir()
repo_root <- dirname(config_dir)

source(file.path(config_dir, "shiny", "app_panels.R"))
source(file.path(config_dir, "shiny", "app_server.R"))

# Define UI for application
ui <- fluidPage(
  # Application title
  theme = shinythemes::shinytheme("simplex"),
  titlePanel("Petrisnake config: Generate configuration yaml file"),
  navlistPanel(id = "mainNav", widths = c(2, 10),
               general_option_panel,
               sample_configuration_panel,
               load_yaml_panel,
               save_yaml_panel
               )
)

# Run the application
shinyApp(ui = ui, server = server)
