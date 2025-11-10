# Define server logic
server <- function(input, output, session) {
    sample_counter <- reactiveVal(0L)
    active_tabs <-
        prefix <- reactive(ifelse(input$prefix |> endsWith("/") || input$prefix == "", input$prefix, paste0(input$prefix, "/")))
    output$path_reference_genome <- renderText(glue("Relative path: {prefix()}{input$reference_genome}"))
    output$path_reference_annotation <- renderText(glue("Relative path: {prefix()}{input$reference_annotation}"))

    observe(
        #Creates new tab
        {
            this_sample_number <- sample_counter() + 1L

            override_params <- c("bc_cutoff", "chunk_size", "feature_tag", "gene_attribute") |> create_self_naming_list()
            override_elements <- c("specified", "field", "value") |> create_self_naming_list()

            basic_ids <- c("tab", "name", "title", "prefix", "suffix", "forward_suffix", "reverse_suffix") |>
              create_self_naming_list()  |>
              map(\(id) "sample_{id}_{this_sample_number}" |> glue())

            override_ids <- map(
                override_params,
                \(parameter) map(
                    override_elements,
                    \(element) glue("sample_{parameter}_{element}_{this_sample_number}")
                )
            )

            ids <- c(basic_ids, override_ids)

            appendTab(
                session = session,
                inputId = "samples_panel",
                tabPanel(
                    title = uiOutput(ids$title),
                    value = ids$tab,
                    wellPanel(
                        h4("Path options"),
                        textInput(
                            inputId = ids$name,
                            label = "Sample name",
                            value = glue("Sample {this_sample_number}")
                        ),
                        textInput(
                            inputId = ids$prefix,
                            label = "Common prefix for all read files in sample",
                            value = ""
                        ),
                        textInput(
                            inputId = ids$suffix,
                            label = "Common suffix for all read files in sample",
                            value = ""
                        ),
                        textInput(
                            inputId = ids$forward_suffix,
                            label = "Suffix for forward read files",
                            value = ""
                        ),
                        textInput(
                            inputId = ids$reverse_suffix,
                            label = "Suffix for reverse read files",
                            value = ""
                        )
                    ),
                    wellPanel(
                        h4("Parameter options"),
                        checkboxInput(ids$bc_cutoff$specified, "Use same barcode cutoff as in the general settings?", value = TRUE),
                        uiOutput(ids$bc_cutoff$field),
                        checkboxInput(ids$chunk_size$specified, "Use same streaming chunk size as in the general settings?", value = TRUE),
                        uiOutput(ids$chunk_size$field),
                        checkboxInput(ids$feature_tag$specified, "Use same GTF feature tag as in the general settings?", value = TRUE),
                        uiOutput(ids$feature_tag$field),
                        checkboxInput(ids$gene_attribute$specified, "Use same GTF gene attribute tag as in the general settings?", value = TRUE),
                        uiOutput(ids$gene_attribute$field)
                    )
                )
            )

            output[[ids$title]] <- renderText({
                req(input[[ids$name]])
                input[[ids$name]]
            })

            output[[ids$bc_cutoff$field]] <- renderUI({
                req(input$bc_cutoff)
                if (isTRUE(input[[ids$bc_cutoff$specified]])) {
                    glue("Keeping {input$bc_cutoff} barcodes") |> h6()
                } else {
                    numericInput(ids$bc_cutoff$value, "Number of barcode combinations to use:", value = 1e4, step = 1000)
                }
            })

            output[[ids$chunk_size$field]] <- renderUI({
                req(input$chunk_size)
                if (isTRUE(input[[ids$chunk_size$specified]])) {
                    glue("Chunk size for streaming: {input$chunk_size} reads") |> h6()
                } else {
                    numericInput(ids$chunk_size$value, "Streaming chunk size:", value = 2e5, step = 1e4)
                }
            })

            output[[ids$feature_tag$field]] <- renderUI({
                req(input$feature_tag)
                if (isTRUE(input[[ids$feature_tag$specified]])) {
                    glue("GTF file feature tag: {input$feature_tag}") |> h6()
                } else {
                    textInput(ids$feature_tag$value, "GTF file feature tag", value = "Coding_or_RNA")
                }
            })

            output[[ids$gene_attribute$field]] <- renderUI({
                req(input$gene_id_attribute)
                if (isTRUE(input[[ids$gene_attribute$specified]])) {
                    glue("GTF file feature tag: {input$gene_id_attribute}") |> h6()
                } else {
                    textInput(ids$gene_attribute$value, "GTF file tag for gene identifiers", value = "name")
                }
            })

            sample_counter(this_sample_number)
        }
    ) |> bindEvent(input$add_sample)

    observe({
        req(input$samples_panel)
        removeTab(
            session = session,
            inputId = "samples_panel",
            target = input$samples_panel
        )
    }) |> bindEvent(input$remove_sample)

    
    observeEvent(input$SaveYAML, {
        write_yaml(
            x = makeYAML(input),
            file = input$savePath
        )
    })

    output$downloadData <- downloadHandler(
        filename = function() {
            paste(input$runID, ".yaml", sep = "")
        },
        content = function(file) {
            write_yaml(
                makeYAML(input),
                file
            )
            # write.csv(datasetInput(), file, row.names = FALSE)
        }
    )

    makeYAML <- function(input) {
        # collect fastq file information
        seqf <- list()
        for (i in 1:input$nfiles) {
            bc_struc <- c(input[[paste0("cDNA_", i)]], input[[paste0("BC_", i)]], input[[paste0("UMI_", i)]])
            names(bc_struc) <- c("cDNA", "BC", "UMI")
            bc_struc <- bc_struc[which(bc_struc != "")]


            if (input$patternsearch == F & input$frameshiftsearch == F) {
                seqf[[i]] <- list(
                    "name" = input[[paste0("fqpath_", i)]],
                    "base_definition" = paste0(names(bc_struc), "(", bc_struc, ")")
                )
            } else {
                if (input$patternsearch == T & substr(input$patternRead, 6, 6) == i) {
                    seqf[[i]] <- list(
                        "name" = input[[paste0("fqpath_", i)]],
                        "base_definition" = paste0(names(bc_struc), "(", bc_struc, ")"),
                        "find_pattern" = input$pattern
                    )
                } else {
                    if (input$frameshiftsearch == T & substr(input$patternRead, 6, 6) == i) {
                        seqf[[i]] <- list(
                            "name" = input[[paste0("fqpath_", i)]],
                            "base_definition" = paste0(names(bc_struc), "(", bc_struc, ")"),
                            "correct_frameshift" = input$pattern
                        )
                    } else {
                        seqf[[i]] <- list(
                            "name" = input[[paste0("fqpath_", i)]],
                            "base_definition" = paste0(names(bc_struc), "(", bc_struc, ")")
                        )
                    }
                }
            }
        }
        names(seqf) <- paste0("file", 1:input$nfiles)

        # collect potential additional fasta files
        if (input$NUMadditionalFA > 0) {
            fa_list <- c()
            for (i in 1:input$NUMadditionalFA) {
                fa_list <- c(fa_list, input[[paste0("FA_", i)]])
            }
        } else {
            fa_list <- NULL
        }

        # decide on barcode param
        # if(input$barcodeChoice)

        # make yaml list
        y <- list(
            "project" = input$runID,
            "sequence_files" = seqf,
            "reference" = list(
                "STAR_index" = input$STARpath,
                "GTF_file" = input$GTFpath,
                "additional_STAR_params" = input$STARparams,
                "additional_files" = fa_list
            ),
            "out_dir" = input$outDir,
            "num_threads" = input$numThreads,
            "mem_limit" = input$memLimit,
            "filter_cutoffs" = list(
                "BC_filter" = list(
                    "num_bases" = input$BCbases,
                    "phred" = input$BCphred
                ),
                "UMI_filter" = list(
                    "num_bases" = input$UMIbases,
                    "phred" = input$UMIphred
                )
            ),
            "barcodes" = list(
                "barcode_num" = input$BCnum,
                "barcode_file" = input$BCfile,
                "barcode_sharing" = input$sharedBC,
                "automatic" = ifelse(input$barcodeChoice == "Automatic", TRUE, FALSE),
                "BarcodeBinning" = input$HamBC,
                "nReadsperCell" = input$nReadsBC,
                "demultiplex" = input$demux
            ),
            "counting_opts" = list(
                "introns" = input$countIntrons,
                "downsampling" = input$downsamp,
                "strand" = as.integer(input$strand),
                "Ham_Dist" = input$HamDist,
                "write_ham" = input$writeHam,
                "velocyto" = input$doVelocity,
                "primaryHit" = input$countPrimary,
                "twoPass" = input$twoPass
            ),
            "make_stats" = input$makeStats,
            "which_Stage" = input$whichStage,
            "Rscript_exec" = input$r_exec,
            "STAR_exec" = input$star_exec,
            "pigz_exec" = input$pigz_exec,
            "samtools_exec" = input$samtools_exec
        )
        return(y)
    }

    observeEvent(input$LoadYAML, {
        if (file.exists(input$inYAML)) {
            print(paste("Loading", input$inYAML))
            loading_func(input$inYAML)
            updateNavlistPanel(session = session, inputId = "mainNav", selected = "Mandatory Parameters")
        } else {
            if (is.null(input$file1)) {
                print("File doesn't exist!")
            } else {
                print(paste("Loading", input$file1$datapath))
                loading_func(input$file1$datapath)
                updateNavlistPanel(session = session, inputId = "mainNav", selected = "Mandatory Parameters")
            }
        }
    })

    loading_func <- function(myYAML) {
        ya <- read_yaml(file = myYAML)
        updateTextInput(session = session, inputId = "runID", value = ya$project)
        updateTextInput(session = session, inputId = "outDir", value = ya$out_dir)
        updateTextInput(session = session, inputId = "STARpath", value = ya$reference$STAR_index)
        updateTextInput(session = session, inputId = "GTFpath", value = ya$reference$GTF_file)
        updateTextInput(session = session, inputId = "STARparams", value = ya$reference$additional_STAR_params)
        updateNumericInput(session = session, inputId = "NUMadditionalFA", value = length(ya$reference$additional_files))
        if (length(ya$reference$additional_files) > 0) {
            for (i in 1:length(ya$reference$additional_files)) {
                updateTextInput(session = session, inputId = paste0("FA_", i), value = ya$reference$additional_files[i])
            }
        }
        updateSliderInput(session = session, inputId = "nfiles", value = length(ya$sequence_files))
        for (i in 1:length(ya$sequence_files)) {
            updateTextInput(session = session, inputId = paste0("fqpath_", i), value = ya$sequence_files[[i]]$name)

            if (any(grepl("BC", ya$sequence_files[[i]]$base_definition))) {
                bc_string <- grep("BC", ya$sequence_files[[i]]$base_definition, value = T)
                bc_string <- substr(bc_string, start = 4, stop = nchar(bc_string) - 1)
                updateTextInput(session = session, inputId = paste0("BC_", i), value = bc_string)
            }
            if (any(grepl("UMI", ya$sequence_files[[i]]$base_definition))) {
                umi_string <- grep("UMI", ya$sequence_files[[i]]$base_definition, value = T)
                umi_string <- substr(umi_string, start = 5, stop = nchar(umi_string) - 1)
                updateTextInput(session = session, inputId = paste0("UMI_", i), value = umi_string)
            }
            if (any(grepl("cDNA", ya$sequence_files[[i]]$base_definition))) {
                cdna_string <- grep("cDNA", ya$sequence_files[[i]]$base_definition, value = T)
                cdna_string <- substr(cdna_string, start = 6, stop = nchar(cdna_string) - 1)
                updateTextInput(session = session, inputId = paste0("cDNA_", i), value = cdna_string)
            }
            if (length(ya$sequence_files[[i]]$find_pattern) == 1) {
                updateCheckboxInput(session = session, inputId = "patternsearch", value = T)
                updateSelectInput(session = session, inputId = "patternRead", choices = paste("Read", 1:length(ya$sequence_files)), selected = paste("Read", i))
                updateTextInput(session = session, inputId = "pattern", value = ya$sequence_files[[i]]$find_pattern)
            }
            if (length(ya$sequence_files[[i]]$correct_frameshift) == 1) {
                updateCheckboxInput(session = session, inputId = "frameshiftsearch", value = T)
                updateSelectInput(session = session, inputId = "patternRead", choices = paste("Read", 1:length(ya$sequence_files)), selected = paste("Read", i))
                updateTextInput(session = session, inputId = "pattern", value = ya$sequence_files[[i]]$correct_frameshift)
            }
        }

        updateNumericInput(session = session, inputId = "BCbases", value = ya$filter_cutoffs$BC_filter$num_bases)
        updateNumericInput(session = session, inputId = "BCphred", value = ya$filter_cutoffs$BC_filter$phred)
        updateNumericInput(session = session, inputId = "UMIbases", value = ya$filter_cutoffs$UMI_filter$num_bases)
        updateNumericInput(session = session, inputId = "UMIphred", value = ya$filter_cutoffs$UMI_filter$phred)
        updateNumericInput(session = session, inputId = "numThreads", value = ya$num_threads)
        updateNumericInput(session = session, inputId = "memLimit", value = ya$mem_limit)
        updateCheckboxInput(session = session, inputId = "makeStats", value = ya$make_stats)
        updateSelectInput(session = session, inputId = "whichStage", selected = ya$which_Stage)
        updateCheckboxInput(session = session, inputId = "countIntrons", value = ya$counting_opts$introns)
        updateTextInput(session = session, inputId = "downsamp", value = ya$counting_opts$downsampling)
        updateSelectInput(session = session, inputId = "strand", selected = ya$counting_opts$strand)
        updateNumericInput(session = session, inputId = "HamDist", value = ya$counting_opts$Ham_Dist)
        updateCheckboxInput(session = session, inputId = "writeHam", value = ya$counting_opts$write_ham)
        updateCheckboxInput(session = session, inputId = "doVelocity", value = ya$counting_opts$velocyto)
        updateCheckboxInput(session = session, inputId = "countPrimary", value = ya$counting_opts$primaryHit)
        updateCheckboxInput(session = session, inputId = "twoPass", value = ya$counting_opts$twoPass)
        if (is.null(ya$barcodes$barcode_num) & ya$barcodes$automatic == TRUE) {
            updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Automatic")
            updateTextInput(session = session, inputId = "BCfile", value = ya$barcodes$barcode_file)
        }
        if (!is.null(ya$barcodes$barcode_file) & ya$barcodes$automatic == FALSE) {
            updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Barcode whitelist")
            updateTextInput(session = session, inputId = "BCfile", value = ya$barcodes$barcode_file)
        }
        if (!is.null(ya$barcodes$barcode_num)) {
            updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Number of top Barcodes")
            updateNumericInput(session = session, inputId = "BCnum", value = ya$barcodes$barcode_num)
        }

        updateNumericInput(session = session, inputId = "HamBC", value = ya$barcodes$BarcodeBinning)
        updateNumericInput(session = session, inputId = "nReadsBC", value = ya$barcodes$nReadsperCell)
        updateTextInput(session = session, inputId = "sharedBC", value = ya$barcodes$barcode_sharing)
        updateCheckboxInput(session = session, inputId = "demux", value = ya$barcodes$demultiplex)

        if (!is.null(ya$read_layout)) {
            updateSelectInput(session = session, inputId = "layout", selected = ya$read_layout)
        }

        updateTextInput(session = session, inputId = "r_exec", value = ya$Rscript_exec)
        updateTextInput(session = session, inputId = "samtools_exec", value = ya$samtools_exec)
        updateTextInput(session = session, inputId = "pigz_exec", value = ya$pigz_exec)
        updateTextInput(session = session, inputId = "star_exec", value = ya$STAR_exec)
    }
}

create_self_naming_list <- function(x) {
    res <- as.list(x)
    names(res) <- x
    res
}
