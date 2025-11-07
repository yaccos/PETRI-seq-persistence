library(shiny)
library(DT)
library(yaml)

config_path <- file.path("config.yaml")

read_config <- function(path) {
  if (!file.exists(path)) {
    list(
      prefix = "data/",
      samples = list(),
      reference_genome = "",
      reference_annotation = "",
      feature_tag = "",
      gene_id_attribute = "",
      streaming_chunk_size = 4000L
    )
  } else {
    yaml.load_file(path)
  }
}

samples_to_df <- function(samples_list) {
  if (!length(samples_list)) return(data.frame(
    sample = character(), bc_cutoff = integer(), prefix = character(),
    lanes = character(), forward_suffix = character(),
    reverse_suffix = character(), suffix = character(), stringsAsFactors = FALSE
  ))
  do.call(rbind, lapply(names(samples_list), function(nm) {
    smp <- samples_list[[nm]]
    data.frame(
      sample = nm,
      bc_cutoff = smp$bc_cutoff,
      prefix = smp$prefix,
      lanes = paste(names(smp$lanes), smp$lanes, sep = ":", collapse = ","),
      forward_suffix = smp$forward_suffix,
      reverse_suffix = smp$reverse_suffix,
      suffix = smp$suffix,
      stringsAsFactors = FALSE
    )
  }))
}

df_to_samples <- function(df) {
  out <- vector("list", nrow(df))
  names(out) <- df$sample
  for (i in seq_len(nrow(df))) {
    lane_pairs <- strsplit(df$lanes[i], ",")[[1]]
    lane_pairs <- lane_pairs[nzchar(lane_pairs)]
    lanes <- setNames(
      sub("^[^:]+:", "", lane_pairs),
      sub(":.*$", "", lane_pairs)
    )
    out[[i]] <- list(
      bc_cutoff = as.integer(df$bc_cutoff[i]),
      prefix = df$prefix[i],
      lanes = as.list(lanes),
      forward_suffix = df$forward_suffix[i],
      reverse_suffix = df$reverse_suffix[i],
      suffix = df$suffix[i]
    )
  }
  out
}

save_config <- function(cfg, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(as.yaml(cfg), path)
}

ui <- fluidPage(
  titlePanel("PETRI-seq Configurator"),
  sidebarLayout(
    sidebarPanel(
      textInput("prefix", "Data prefix", ""),
      textInput("ref_genome", "Reference genome", ""),
      textInput("ref_annotation", "Reference annotation", ""),
      textInput("feature_tag", "Feature tag", ""),
      textInput("gene_id_attr", "Gene ID attribute", ""),
      numericInput("chunk_size", "Streaming chunk size", value = 4000, min = 1, step = 100),
      hr(),
      actionButton("add_sample", "Add sample"),
      actionButton("remove_sample", "Remove selected"),
      hr(),
      actionButton("save_cfg", "Save config"),
      hr(),
    ),
    mainPanel(
      DTOutput("samples_table"),
      hr(),
      verbatimTextOutput("status")
    )
  )
)

server <- function(input, output, session) {
  cfg <- reactiveVal(read_config(config_path))
  samples_df <- reactiveVal(samples_to_df(cfg()$samples))
  status_txt <- reactiveVal("Ready.")

  observe({
    current <- cfg()
    updateTextInput(session, "prefix", value = current$prefix)
    updateTextInput(session, "ref_genome", value = current$reference_genome)
    updateTextInput(session, "ref_annotation", value = current$reference_annotation)
    updateTextInput(session, "feature_tag", value = current$feature_tag)
    updateTextInput(session, "gene_id_attr", value = current$gene_id_attribute)
    updateNumericInput(session, "chunk_size", value = current$streaming_chunk_size)
  })

  output$samples_table <- renderDT({
    datatable(
      samples_df(),
      editable = TRUE,
      selection = "single",
      options = list(pageLength = 10, dom = "tp")
    )
  }, server = FALSE)

  observeEvent(input$add_sample, {
    df <- samples_df()
    df <- rbind(df, data.frame(
      sample = paste0("sample_", nrow(df) + 1),
      bc_cutoff = 1000L,
      prefix = "",
      lanes = "L1:L001",
      forward_suffix = "_R1",
      reverse_suffix = "_R2",
      suffix = "_001.fastq.gz",
      stringsAsFactors = FALSE
    ))
    samples_df(df)
  })

  observeEvent(input$remove_sample, {
    idx <- input$samples_table_rows_selected
    if (length(idx)) {
      df <- samples_df()
      df <- df[-idx, ]
      samples_df(df)
    }
  })

  observeEvent(input$samples_table_cell_edit, {
    info <- input$samples_table_cell_edit
    df <- samples_df()
    df[info$row, info$col] <- info$value
    samples_df(df)
  })

  observeEvent(input$save_cfg, {
    new_cfg <- list(
      prefix = input$prefix,
      samples = df_to_samples(samples_df()),
      reference_genome = input$ref_genome,
      reference_annotation = input$ref_annotation,
      feature_tag = input$feature_tag,
      gene_id_attribute = input$gene_id_attr,
      streaming_chunk_size = as.integer(input$chunk_size)
    )
    save_config(new_cfg, config_path)
    cfg(new_cfg)
    status_txt("Configuration saved.")
  })

  output$status <- renderText(status_txt())
}

shinyApp(ui, server)

