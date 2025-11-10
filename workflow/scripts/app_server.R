# Define server logic
server <- function(input, output, session) {
    sample_counter <- reactiveVal(0L)
    prefix <- reactive(ifelse(input$prefix |> endsWith("/") || input$prefix == "", input$prefix, paste0(input$prefix, "/")))

    sample_registry <- reactiveVal(list())

    register_sample <- function(ids) {
        entries <- sample_registry()
        entries[[ids$tab]] <- ids
        sample_registry(entries)
    }

    unregister_sample <- function(tab_id) {
        entries <- sample_registry()
        if (!is.null(entries[[tab_id]])) {
            entries[[tab_id]] <- NULL
            sample_registry(entries)
        }
    }

    ensure_unique_name <- function(name, used_names) {
        candidate <- name
        suffix <- 2L
        while (candidate %in% used_names) {
            candidate <- paste0(name, "_", suffix)
            suffix <- suffix + 1L
        }
        candidate
    }

    `%||%` <- function(x, y) {
        if (is.null(x)) y else x
    }

    run_after_flush <- function(callback) {
        session$onFlushed(callback, once = TRUE)
    }

    coerce_numeric <- function(value) {
        if (is.null(value)) {
            return(NULL)
        }
        if (length(value) == 0) {
            return(NULL)
        }
        if (is.numeric(value) && length(value) == 1L && !is.na(value)) {
            return(value)
        }
        numeric_value <- suppressWarnings(as.numeric(value))
        if (length(numeric_value) == 1L && !is.na(numeric_value)) {
            numeric_value
        } else {
            NULL
        }
    }

    input_value <- function(id, default = "") {
        value <- input[[id]]
        if (is.null(value)) {
            default
        } else {
            value
        }
    }
    output$path_reference_genome <- renderText(glue("Relative path: {prefix()}{input$reference_genome}"))
    output$path_reference_annotation <- renderText(glue("Relative path: {prefix()}{input$reference_annotation}"))

    add_sample_tab <- function(initial = NULL) {
        this_sample_number <- sample_counter() + 1L

        initial <- initial %||% list()
        loading_from_config <- isTRUE(initial$.__from_config__)
        initial$.__from_config__ <- NULL

        normalize_text <- function(value) {
            if (is.null(value)) {
                return(NULL)
            }
            trimws(as.character(value))
        }

        sample_label <- normalize_text(initial$name) %||% glue("Sample {this_sample_number}")
        sample_prefix <- normalize_text(initial$prefix) %||% ""

        lane_defaults <- initial$lanes
        if (is.null(lane_defaults)) {
            lane_defaults <- list()
        }
        lane_defaults <- as.list(lane_defaults)
        lane_default_names <- names(lane_defaults)
        if (is.null(lane_default_names)) {
            lane_default_names <- character(length(lane_defaults))
        } else {
            lane_default_names <- trimws(as.character(lane_default_names))
        }
        if (length(lane_defaults) > 0) {
            lane_defaults <- map(lane_defaults, normalize_text)
            names(lane_defaults) <- lane_default_names
        }
        lane_count_default <- max(1L, length(lane_defaults))

        override_lookup <- list(
            bc_cutoff = "bc_cutoff",
            chunk_size = "streaming_chunk_size",
            feature_tag = "feature_tag",
            gene_attribute = "gene_id_attribute",
            suffix = "suffix",
            forward_suffix = "forward_suffix",
            reverse_suffix = "reverse_suffix",
            reference_genome = "reference_genome",
            reference_annotation = "reference_annotation"
        )

        override_defaults <- imap(override_lookup, function(source_name, key) {
            raw_value <- initial[[source_name]]
            if (is.null(raw_value)) {
                return(list(use_general = TRUE, value = NULL))
            }
            if (key %in% c("bc_cutoff", "chunk_size")) {
                normalized <- coerce_numeric(raw_value)
            } else {
                normalized <- trimws(as.character(raw_value))
            }
            list(use_general = FALSE, value = normalized)
        })

        override_params <- c(
            "bc_cutoff",
            "chunk_size",
            "feature_tag",
            "gene_attribute",
            "suffix",
            "forward_suffix",
            "reverse_suffix",
            "reference_genome",
            "reference_annotation"
        ) |> create_self_naming_list()
        override_elements <- c("specified", "field", "value") |> create_self_naming_list()

        basic_ids <- c("tab", "name", "title", "prefix", "lane_count", "lanes") |>
            create_self_naming_list() |>
            map(\(id) "sample_{id}_{this_sample_number}" |> glue())

        override_ids <- map(
            override_params,
            \(parameter) map(
                override_elements,
                \(element) glue("sample_{parameter}_{element}_{this_sample_number}")
            )
        )

        ids <- c(basic_ids, override_ids)

        resolve_override <- function(key, general_id, default = "") {
            use_general <- input[[ids[[key]]$specified]]
            if (is.null(use_general)) {
                use_general <- override_defaults[[key]]$use_general
            }
            if (isTRUE(use_general)) {
                input_value(general_id, default)
            } else {
                input_value(ids[[key]]$value, default)
            }
        }

        appendTab(
            session = session,
            inputId = "samples_panel",
            tabPanel(
                title = uiOutput(ids$title),
                value = ids$tab,
                fluidRow(
                    column(
                        width = 6,
                        wellPanel(
                            h4("Path options"),
                            textInput(
                                inputId = ids$name,
                                label = "Sample name",
                                value = sample_label
                            ),
                            textInput(
                                inputId = ids$prefix,
                                label = "Common prefix for all read files in sample",
                                value = sample_prefix
                            ),
                            checkboxInput(ids$suffix$specified, "Use common suffix from general settings?", value = override_defaults$suffix$use_general),
                            uiOutput(ids$suffix$field),
                            checkboxInput(ids$forward_suffix$specified, "Use forward read suffix from general settings?", value = override_defaults$forward_suffix$use_general),
                            uiOutput(ids$forward_suffix$field),
                            checkboxInput(ids$reverse_suffix$specified, "Use reverse read suffix from general settings?", value = override_defaults$reverse_suffix$use_general),
                            uiOutput(ids$reverse_suffix$field)
                        ),
                        wellPanel(
                            h4("Parameter options"),
                            checkboxInput(ids$bc_cutoff$specified, "Use same barcode cutoff as in the general settings?", value = override_defaults$bc_cutoff$use_general),
                            uiOutput(ids$bc_cutoff$field),
                            checkboxInput(ids$chunk_size$specified, "Use same streaming chunk size as in the general settings?", value = override_defaults$chunk_size$use_general),
                            uiOutput(ids$chunk_size$field),
                            checkboxInput(ids$feature_tag$specified, "Use same GTF feature tag as in the general settings?", value = override_defaults$feature_tag$use_general),
                            uiOutput(ids$feature_tag$field),
                            checkboxInput(ids$gene_attribute$specified, "Use same GTF gene attribute tag as in the general settings?", value = override_defaults$gene_attribute$use_general),
                            uiOutput(ids$gene_attribute$field)
                        ),
                        wellPanel(
                            h4("Reference files"),
                            checkboxInput(ids$reference_genome$specified, "Use reference genome from general settings?", value = override_defaults$reference_genome$use_general),
                            uiOutput(ids$reference_genome$field),
                            checkboxInput(ids$reference_annotation$specified, "Use reference annotation from general settings?", value = override_defaults$reference_annotation$use_general),
                            uiOutput(ids$reference_annotation$field)
                        )
                    ),
                    column(
                        width = 6,
                        wellPanel(
                            h4("Lane configuration"),
                            numericInput(ids$lane_count, "Number of lanes", value = lane_count_default, min = 1, step = 1),
                            uiOutput(ids$lanes)
                        )
                    )
                )
            )
        )

        output[[ids$title]] <- renderText({
            req(input[[ids$name]])
            input[[ids$name]]
        })

        output[[ids$bc_cutoff$field]] <- renderUI({
            req(input$bc_cutoff)
            specified <- input[[ids$bc_cutoff$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$bc_cutoff$use_general
            }
            if (isTRUE(specified)) {
                glue("Keeping {input$bc_cutoff} barcodes") |> h6()
            } else {
                numericInput(
                    ids$bc_cutoff$value,
                    "Number of barcode combinations to use:",
                    value = override_defaults$bc_cutoff$value %||% input$bc_cutoff,
                    step = 1000
                )
            }
        })

        output[[ids$chunk_size$field]] <- renderUI({
            req(input$chunk_size)
            specified <- input[[ids$chunk_size$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$chunk_size$use_general
            }
            if (isTRUE(specified)) {
                glue("Chunk size for streaming: {input$chunk_size} reads") |> h6()
            } else {
                numericInput(
                    ids$chunk_size$value,
                    "Streaming chunk size:",
                    value = override_defaults$chunk_size$value %||% input$chunk_size,
                    step = 1e4
                )
            }
        })

        output[[ids$feature_tag$field]] <- renderUI({
            req(input$feature_tag)
            specified <- input[[ids$feature_tag$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$feature_tag$use_general
            }
            if (isTRUE(specified)) {
                glue("GTF file feature tag: {input$feature_tag}") |> h6()
            } else {
                textInput(
                    ids$feature_tag$value,
                    "GTF file feature tag",
                    value = override_defaults$feature_tag$value %||% input_value("feature_tag", "")
                )
            }
        })

        output[[ids$gene_attribute$field]] <- renderUI({
            req(input$gene_id_attribute)
            specified <- input[[ids$gene_attribute$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$gene_attribute$use_general
            }
            if (isTRUE(specified)) {
                glue("GTF file feature tag: {input$gene_id_attribute}") |> h6()
            } else {
                textInput(
                    ids$gene_attribute$value,
                    "GTF file tag for gene identifiers",
                    value = override_defaults$gene_attribute$value %||% input_value("gene_id_attribute", "")
                )
            }
        })

        output[[ids$suffix$field]] <- renderUI({
            suffix_value <- input_value("suffix")
            specified <- input[[ids$suffix$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$suffix$use_general
            }
            if (isTRUE(specified)) {
                glue("Common suffix: {suffix_value}") |> h6()
            } else {
                textInput(
                    ids$suffix$value,
                    "Common suffix for all read files",
                    value = override_defaults$suffix$value %||% suffix_value
                )
            }
        })

        output[[ids$forward_suffix$field]] <- renderUI({
            suffix_value <- input_value("forward_suffix")
            specified <- input[[ids$forward_suffix$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$forward_suffix$use_general
            }
            if (isTRUE(specified)) {
                glue("Forward read suffix: {suffix_value}") |> h6()
            } else {
                textInput(
                    ids$forward_suffix$value,
                    "Suffix for forward read files",
                    value = override_defaults$forward_suffix$value %||% suffix_value
                )
            }
        })

        output[[ids$reverse_suffix$field]] <- renderUI({
            suffix_value <- input_value("reverse_suffix")
            specified <- input[[ids$reverse_suffix$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$reverse_suffix$use_general
            }
            if (isTRUE(specified)) {
                glue("Reverse read suffix: {suffix_value}") |> h6()
            } else {
                textInput(
                    ids$reverse_suffix$value,
                    "Suffix for reverse read files",
                    value = override_defaults$reverse_suffix$value %||% suffix_value
                )
            }
        })

        reference_genome_path_id <- glue("{ids$reference_genome$value}_path")

        output[[reference_genome_path_id]] <- renderText({
            req(!isTRUE(input[[ids$reference_genome$specified]]))
            sample_genome_value <- input_value(ids$reference_genome$value, "")
            glue("Relative path: {prefix()}{sample_genome_value}")
        })

        output[[ids$reference_genome$field]] <- renderUI({
            genome_value <- input_value("reference_genome")
            specified <- input[[ids$reference_genome$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$reference_genome$use_general
            }
            if (isTRUE(specified)) {
                if (nzchar(genome_value)) {
                    glue("Reference genome: {prefix()}{genome_value}") |> h6()
                } else {
                    h6("Reference genome not set")
                }
            } else {
                existing_value <- isolate(input[[ids$reference_genome$value]])
                sample_genome_value <- if (is.null(existing_value)) {
                    override_defaults$reference_genome$value %||% genome_value
                } else {
                    existing_value
                }
                tagList(
                    textInput(ids$reference_genome$value, "Reference genome for this sample", value = sample_genome_value),
                    h6(textOutput(reference_genome_path_id))
                )
            }
        })

        reference_annotation_path_id <- glue("{ids$reference_annotation$value}_path")

        output[[reference_annotation_path_id]] <- renderText({
            req(!isTRUE(input[[ids$reference_annotation$specified]]))
            sample_annotation_value <- input_value(ids$reference_annotation$value, "")
            glue("Relative path: {prefix()}{sample_annotation_value}")
        })

        output[[ids$reference_annotation$field]] <- renderUI({
            annotation_value <- input_value("reference_annotation")
            specified <- input[[ids$reference_annotation$specified]]
            if (is.null(specified)) {
                specified <- override_defaults$reference_annotation$use_general
            }
            if (isTRUE(specified)) {
                if (nzchar(annotation_value)) {
                    glue("Reference annotation: {prefix()}{annotation_value}") |> h6()
                } else {
                    h6("Reference annotation not set")
                }
            } else {
                existing_value <- isolate(input[[ids$reference_annotation$value]])
                sample_annotation_value <- if (is.null(existing_value)) {
                    override_defaults$reference_annotation$value %||% annotation_value
                } else {
                    existing_value
                }
                tagList(
                    textInput(ids$reference_annotation$value, "Reference annotation for this sample", value = sample_annotation_value),
                    h6(textOutput(reference_annotation_path_id))
                )
            }
        })

        output[[ids$lanes]] <- renderUI({
            lane_total <- suppressWarnings(as.integer(input_value(ids$lane_count, lane_count_default)))
            if (is.na(lane_total) || lane_total < 1) {
                return(h6("No lanes configured"))
            }

            lane_ui <- lapply(seq_len(lane_total), function(idx) {
                lane_name_id <- glue("{ids$tab}_lane_name_{idx}")
                lane_identifier_id <- glue("{ids$tab}_lane_identifier_{idx}")
                forward_path_id <- glue("{lane_identifier_id}_forward_path")
                reverse_path_id <- glue("{lane_identifier_id}_reverse_path")
                lane_name_value <- isolate(input[[lane_name_id]])
                if (is.null(lane_name_value)) {
                    lane_name_value <- lane_default_names[[idx]] %||% glue("L{idx}")
                }
                lane_identifier_value <- isolate(input[[lane_identifier_id]])
                if (is.null(lane_identifier_value)) {
                    lane_identifier_value <- lane_defaults[[idx]] %||% ""
                }

                local({
                    lane_identifier_id_local <- lane_identifier_id
                    forward_path_id_local <- forward_path_id
                    reverse_path_id_local <- reverse_path_id

                    output[[forward_path_id_local]] <- renderText({
                        lane_identifier_current <- input[[lane_identifier_id_local]]
                        if (is.null(lane_identifier_current)) {
                            lane_identifier_current <- ""
                        }
                        suffix_value <- resolve_override("suffix", "suffix")
                        forward_suffix_value <- resolve_override("forward_suffix", "forward_suffix")
                        sample_prefix_value <- input_value(ids$prefix, "")
                        paste0(prefix(), sample_prefix_value, lane_identifier_current, forward_suffix_value, suffix_value)
                    })

                    output[[reverse_path_id_local]] <- renderText({
                        lane_identifier_current <- input[[lane_identifier_id_local]]
                        if (is.null(lane_identifier_current)) {
                            lane_identifier_current <- ""
                        }
                        suffix_value <- resolve_override("suffix", "suffix")
                        reverse_suffix_value <- resolve_override("reverse_suffix", "reverse_suffix")
                        sample_prefix_value <- input_value(ids$prefix, "")
                        paste0(prefix(), sample_prefix_value, lane_identifier_current, reverse_suffix_value, suffix_value)
                    })
                })

                div(
                    class = "lane-entry",
                    textInput(lane_name_id, glue("Lane {idx} name"), value = lane_name_value),
                    textInput(lane_identifier_id, glue("Lane {idx} identifier"), value = lane_identifier_value),
                    tags$p("Forward read path ", textOutput(forward_path_id, inline = TRUE, container = tags$code)),
                    tags$p("Reverse read path ", textOutput(reverse_path_id, inline = TRUE, container = tags$code))
                )
            })

            tagList(lane_ui)
        })

        register_sample(ids)
        sample_counter(this_sample_number)

        if (loading_from_config) {
            run_after_flush(local({
                ids_local <- ids
                lane_defaults_local <- lane_defaults
                lane_default_names_local <- lane_default_names
                lane_count_local <- lane_count_default
                override_defaults_local <- override_defaults
                sample_label_local <- sample_label
                sample_prefix_local <- sample_prefix
                function() {
                    updateTextInput(session, ids_local$name, value = sample_label_local)
                    updateTextInput(session, ids_local$prefix, value = sample_prefix_local)
                    updateNumericInput(session, ids_local$lane_count, value = lane_count_local)
                    if (lane_count_local > 0) {
                        for (idx in seq_len(lane_count_local)) {
                            lane_name_id <- glue("{ids_local$tab}_lane_name_{idx}")
                            lane_identifier_id <- glue("{ids_local$tab}_lane_identifier_{idx}")
                            lane_name_value <- lane_default_names_local[[idx]] %||% glue("L{idx}")
                            lane_identifier_value <- lane_defaults_local[[idx]] %||% ""
                            updateTextInput(session, lane_name_id, value = lane_name_value)
                            updateTextInput(session, lane_identifier_id, value = lane_identifier_value)
                        }
                    }

                    updateCheckboxInput(session, ids_local$bc_cutoff$specified, value = override_defaults_local$bc_cutoff$use_general)
                    if (!override_defaults_local$bc_cutoff$use_general && !is.null(override_defaults_local$bc_cutoff$value)) {
                        updateNumericInput(session, ids_local$bc_cutoff$value, value = override_defaults_local$bc_cutoff$value)
                    }

                    updateCheckboxInput(session, ids_local$chunk_size$specified, value = override_defaults_local$chunk_size$use_general)
                    if (!override_defaults_local$chunk_size$use_general && !is.null(override_defaults_local$chunk_size$value)) {
                        updateNumericInput(session, ids_local$chunk_size$value, value = override_defaults_local$chunk_size$value)
                    }

                    updateCheckboxInput(session, ids_local$feature_tag$specified, value = override_defaults_local$feature_tag$use_general)
                    if (!override_defaults_local$feature_tag$use_general && !is.null(override_defaults_local$feature_tag$value)) {
                        updateTextInput(session, ids_local$feature_tag$value, value = override_defaults_local$feature_tag$value)
                    }

                    updateCheckboxInput(session, ids_local$gene_attribute$specified, value = override_defaults_local$gene_attribute$use_general)
                    if (!override_defaults_local$gene_attribute$use_general && !is.null(override_defaults_local$gene_attribute$value)) {
                        updateTextInput(session, ids_local$gene_attribute$value, value = override_defaults_local$gene_attribute$value)
                    }

                    updateCheckboxInput(session, ids_local$suffix$specified, value = override_defaults_local$suffix$use_general)
                    if (!override_defaults_local$suffix$use_general && !is.null(override_defaults_local$suffix$value)) {
                        updateTextInput(session, ids_local$suffix$value, value = override_defaults_local$suffix$value)
                    }

                    updateCheckboxInput(session, ids_local$forward_suffix$specified, value = override_defaults_local$forward_suffix$use_general)
                    if (!override_defaults_local$forward_suffix$use_general && !is.null(override_defaults_local$forward_suffix$value)) {
                        updateTextInput(session, ids_local$forward_suffix$value, value = override_defaults_local$forward_suffix$value)
                    }

                    updateCheckboxInput(session, ids_local$reverse_suffix$specified, value = override_defaults_local$reverse_suffix$use_general)
                    if (!override_defaults_local$reverse_suffix$use_general && !is.null(override_defaults_local$reverse_suffix$value)) {
                        updateTextInput(session, ids_local$reverse_suffix$value, value = override_defaults_local$reverse_suffix$value)
                    }

                    updateCheckboxInput(session, ids_local$reference_genome$specified, value = override_defaults_local$reference_genome$use_general)
                    if (!override_defaults_local$reference_genome$use_general && !is.null(override_defaults_local$reference_genome$value)) {
                        updateTextInput(session, ids_local$reference_genome$value, value = override_defaults_local$reference_genome$value)
                    }

                    updateCheckboxInput(session, ids_local$reference_annotation$specified, value = override_defaults_local$reference_annotation$use_general)
                    if (!override_defaults_local$reference_annotation$use_general && !is.null(override_defaults_local$reference_annotation$value)) {
                        updateTextInput(session, ids_local$reference_annotation$value, value = override_defaults_local$reference_annotation$value)
                    }
                }
            }))
        }

        invisible(ids)
    }

    observe({
        add_sample_tab()
    }) |> bindEvent(input$add_sample)

    observe({
        req(input$samples_panel)
        target_tab <- input$samples_panel
        removeTab(
            session = session,
            inputId = "samples_panel",
            target = target_tab
        )
        unregister_sample(target_tab)
    }) |> bindEvent(input$remove_sample)

    collect_samples <- function() {
        entries <- sample_registry()
        if (length(entries) == 0) {
            return(list())
        }

        samples <- list()
        used_sample_names <- character()

        for (tab_id in names(entries)) {
            ids <- entries[[tab_id]]
            raw_name <- trimws(input_value(ids$name, tab_id))
            if (!nzchar(raw_name)) {
                raw_name <- tab_id
            }
            sample_name <- ensure_unique_name(raw_name, used_sample_names)
            used_sample_names <- c(used_sample_names, sample_name)

            sample_prefix <- trimws(input_value(ids$prefix, ""))

            lane_count <- suppressWarnings(as.integer(input_value(ids$lane_count, 0)))
            lane_map <- list()
            lane_names_used <- character()
            if (!is.na(lane_count) && lane_count > 0) {
                for (idx in seq_len(lane_count)) {
                    lane_name_id <- glue("{ids$tab}_lane_name_{idx}")
                    lane_identifier_id <- glue("{ids$tab}_lane_identifier_{idx}")
                    lane_name <- trimws(input_value(lane_name_id, glue("L{idx}")))
                    if (!nzchar(lane_name)) {
                        lane_name <- glue("L{idx}")
                    }
                    lane_name <- ensure_unique_name(lane_name, lane_names_used)
                    lane_names_used <- c(lane_names_used, lane_name)
                    lane_identifier <- trimws(input_value(lane_identifier_id, ""))
                    lane_map[[lane_name]] <- lane_identifier
                }
            }

            sample_entry <- list(
                prefix = sample_prefix,
                lanes = lane_map
            )

            if (!isTRUE(input[[ids$bc_cutoff$specified]])) {
                bc_value <- input[[ids$bc_cutoff$value]]
                if (is.null(bc_value) || is.na(bc_value)) {
                    bc_value <- input$bc_cutoff
                }
                if (!is.null(bc_value) && !is.na(bc_value)) {
                    sample_entry$bc_cutoff <- bc_value
                }
            }

            if (!isTRUE(input[[ids$chunk_size$specified]])) {
                chunk_value <- input[[ids$chunk_size$value]]
                if (is.null(chunk_value) || is.na(chunk_value)) {
                    chunk_value <- input$chunk_size
                }
                if (!is.null(chunk_value) && !is.na(chunk_value)) {
                    sample_entry$streaming_chunk_size <- chunk_value
                }
            }

            if (!isTRUE(input[[ids$feature_tag$specified]])) {
                sample_entry$feature_tag <- trimws(input_value(ids$feature_tag$value, input$feature_tag))
            }

            if (!isTRUE(input[[ids$gene_attribute$specified]])) {
                sample_entry$gene_id_attribute <- trimws(input_value(ids$gene_attribute$value, input$gene_id_attribute))
            }

            if (!isTRUE(input[[ids$suffix$specified]])) {
                sample_entry$suffix <- trimws(input_value(ids$suffix$value, input$suffix))
            }

            if (!isTRUE(input[[ids$forward_suffix$specified]])) {
                sample_entry$forward_suffix <- trimws(input_value(ids$forward_suffix$value, input$forward_suffix))
            }

            if (!isTRUE(input[[ids$reverse_suffix$specified]])) {
                sample_entry$reverse_suffix <- trimws(input_value(ids$reverse_suffix$value, input$reverse_suffix))
            }

            if (!isTRUE(input[[ids$reference_genome$specified]])) {
                sample_entry$reference_genome <- trimws(input_value(ids$reference_genome$value, input$reference_genome))
            }

            if (!isTRUE(input[[ids$reference_annotation$specified]])) {
                sample_entry$reference_annotation <- trimws(input_value(ids$reference_annotation$value, input$reference_annotation))
            }

            samples[[sample_name]] <- sample_entry
        }

        samples
    }

    build_config <- function() {
        config <- list(
            prefix = trimws(input_value("prefix", "")),
            suffix = trimws(input_value("suffix", "")),
            forward_suffix = trimws(input_value("forward_suffix", "")),
            reverse_suffix = trimws(input_value("reverse_suffix", "")),
            reference_genome = trimws(input_value("reference_genome", "")),
            reference_annotation = trimws(input_value("reference_annotation", "")),
            streaming_chunk_size = input$chunk_size,
            bc_cutoff = input$bc_cutoff,
            feature_tag = trimws(input_value("feature_tag", "")),
            gene_id_attribute = trimws(input_value("gene_id_attribute", "")),
            samples = collect_samples()
        )

        if (is.null(config$streaming_chunk_size) || is.na(config$streaming_chunk_size)) {
            config$streaming_chunk_size <- NULL
        }
        if (is.null(config$bc_cutoff) || is.na(config$bc_cutoff)) {
            config$bc_cutoff <- NULL
        }

        Filter(Negate(is.null), config)
    }

    remove_all_samples <- function() {
        current <- sample_registry()
        for (tab_id in names(current)) {
            removeTab(session = session, inputId = "samples_panel", target = tab_id)
        }
        sample_registry(list())
        sample_counter(0L)
    }

    apply_config <- function(cfg) {
        updateTextInput(session, "prefix", value = cfg$prefix %||% "")
        updateTextInput(session, "reference_genome", value = cfg$reference_genome %||% "")
        updateTextInput(session, "reference_annotation", value = cfg$reference_annotation %||% "")
        updateTextInput(session, "suffix", value = cfg$suffix %||% "")
        updateTextInput(session, "forward_suffix", value = cfg$forward_suffix %||% "")
        updateTextInput(session, "reverse_suffix", value = cfg$reverse_suffix %||% "")
        chunk_value <- coerce_numeric(cfg$streaming_chunk_size)
        if (!is.null(chunk_value)) {
            updateNumericInput(session, "chunk_size", value = chunk_value)
        }
        bc_value <- coerce_numeric(cfg$bc_cutoff)
        if (!is.null(bc_value)) {
            updateNumericInput(session, "bc_cutoff", value = bc_value)
        }
        if (!is.null(cfg$feature_tag)) {
            updateTextInput(session, "feature_tag", value = cfg$feature_tag)
        }
        if (!is.null(cfg$gene_id_attribute)) {
            updateTextInput(session, "gene_id_attribute", value = cfg$gene_id_attribute)
        }

        remove_all_samples()

        samples_cfg <- cfg$samples
        if (is.null(names(samples_cfg))) {
            names(samples_cfg) <- paste0("sample_", seq_along(samples_cfg))
        }
        if (length(samples_cfg) == 0) {
            return(invisible())
        }

        for (sample_name in names(samples_cfg)) {
            sample_details <- samples_cfg[[sample_name]]
            if (is.null(sample_details) || !is.list(sample_details)) {
                sample_details <- list()
            }
            sample_details$name <- sample_name
            sample_details$.__from_config__ <- TRUE
            add_sample_tab(initial = sample_details)
        }
    }

    
    observeEvent(input$SaveYAML, {
        req(nzchar(input$savePath))
        config <- build_config()
        tryCatch(
            {
                write_yaml(config, input$savePath)
                showNotification(glue("Saved configuration to {input$savePath}"), type = "message")
            },
            error = function(err) {
                showNotification(glue("Failed to save YAML: {err$message}"), type = "error")
            }
        )
    })

    output$downloadData <- downloadHandler(
        filename = function() {
            if (nzchar(input$savePath)) {
                return(basename(input$savePath))
            }
            paste0("config-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".yaml")
        },
        content = function(file) {
            write_yaml(build_config(), file)
        }
    )

    observeEvent(input$LoadYAML, {
        candidate_path <- NULL
        if (nzchar(input$inYAML) && file.exists(input$inYAML)) {
            candidate_path <- input$inYAML
        } else if (!is.null(input$file1) && file.exists(input$file1$datapath)) {
            candidate_path <- input$file1$datapath
        }

        if (is.null(candidate_path)) {
            showNotification("Please provide a valid YAML path or upload a file.", type = "warning")
            return()
        }

        cfg <- tryCatch(
            read_yaml(candidate_path),
            error = function(err) {
                showNotification(glue("Failed to load YAML: {err$message}"), type = "error")
                NULL
            }
        )
        if (is.null(cfg)) {
            return()
        }

        if (is.null(cfg$samples)) {
            cfg$samples <- list()
        }

        if (!is.list(cfg$samples)) {
            showNotification("The `samples` entry in the YAML file must be a mapping of sample names to settings.", type = "error")
            return()
        }

        apply_config(cfg)
        updateNavlistPanel(session = session, inputId = "mainNav", selected = "General options")
        showNotification(glue("Loaded configuration from {basename(candidate_path)}"), type = "message")
    })
}

create_self_naming_list <- function(x) {
    res <- as.list(x)
    names(res) <- x
    res
}
