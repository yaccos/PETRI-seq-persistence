general_option_panel <- tabPanel(
    "General options",
    fluidRow(
        column(6, wellPanel(
            h4("Data and file input"),
            textInput("prefix", "Directory of input data", value = "data/"),
            textInput("reference_genome", "Reference genome for alignment"),
            textOutput("path_reference_genome") |> h6(),
            # shinyBS::bsTooltip(
            #     id = "reference_genome", title = "Will be overrided for the samples which specify the reference genome differently",
            #     placement = "bottom", trigger = "hover", options = list(container = "body")
            # ),
            textInput("reference_annotation", "Annotation GTF file"),
            textOutput("path_reference_annotation") |> h6(),
            # shinyBS::bsTooltip(
            #     id = "reference_annotation", title = "Will be overrided for the samples which specify the reference annotation differently",
            #     placement = "bottom", trigger = "hover", options = list(container = "body")
            # )
        )),
        column(6, wellPanel(
            h4("Parameter options"),
            numericInput("chunk_size", "Chunk size for data streaming", value = 2e5, step = 10000),
            numericInput("bc_cutoff", "Number of barcode combinations to use", value = 1e4, step = 1000),
            shinyBS::bsTooltip(
                id = "bc_cutoff", title = "Must be specified if you attempting to run the entire pipeline to completion.
                             If you have not decided in advance, run the rule `determine_bc_cutoff` and follow the subsequent instructions.",
                placement = "bottom", trigger = "hover", options = list(container = "body")
            ),
            textInput("feature_tag", "The feature tag in the GTF file to be used for genes", "Coding_or_RNA"),
            textInput("gene_id_attribute", "The attribute in the GTF to be used for the gene identifiers", "name")
        ))
    )
)

sample_configuration_panel <- tabPanel(
    "Sample configuration",
    fluidPage(
        fluidRow(
            actionButton("add_sample", "Add sample", icon = icon("plus-square")),
            actionButton("remove_sample", "Remove selected sample", icon = icon("minus-square"))
      ),
      br(),
      fluidRow(tabsetPanel(
            id = "samples_panel",
        ))
    )
    # # set basic things
    # fluidRow(
    #     column(6, wellPanel(
    #         textInput(inputId = "runID", label = "Name of the run/project:", placeholder = "eg: my_PETRIseq_run")
    #     )),
    #     shinyBS::bsTooltip(
    #         id = "runID", title = "This name will be used to label output files produced by PETRI-seq.",
    #         placement = "bottom", trigger = "hover", options = list(container = "body")
    #     ),
    #     column(6, wellPanel(
    #         textInput(inputId = "outDir", label = "Path to the output directory:", placeholder = "eg: /path/to/output"),
    #         shinyBS::bsTooltip(
    #             id = "outDir", title = "Please remember to give the full path, as relative paths may not work.",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         )
    #     ))
    # ),


    # # slider input for number of fastq reads
    # fluidRow(
    #     h4("Input options:", style = "padding-left: 20px;"),
    #     column(4, wellPanel(
    #         sliderInput("nfiles", "Number of input fastq files:", min = 1, max = 4, value = 2),
    #         shinyBS::bsTooltip(
    #             id = "nfiles", title = "How many reads (including index reads) were obtained by your sequencing layout?",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         ),
    #         checkboxInput("patternsearch", "Search for sequence pattern in reads?", value = F),
    #         checkboxInput("frameshiftsearch", "Correct for frameshift in barcode reads?", value = F),
    #         uiOutput("patternUI"), uiOutput("patternReadUI"),
    #         uiOutput("k"), uiOutput("frameshiftReadUI")
    #     )),
    #     column(8, wellPanel(
    #         uiOutput("fqUI"),
    #         shinyBS::bsTooltip(
    #             id = "fqUI", title = "Please remember to give full paths, as relative paths may not work.",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         )
    #     ))
    # ),
    # fluidRow(
    #     h4("In this section, fill only those fields that fit your input files:"),
    #     column(3, wellPanel(
    #         uiOutput("fqBCui"),
    #         shinyBS::bsTooltip(
    #             id = "fqBCui", title = "List any barcode (BC) ranges to be extracted from the reads. You can also extract several barcode ranges from the same read using comma-separation: 1-6,11-16 ",
    #             placement = "left", trigger = "hover", options = list(container = "body")
    #         )
    #     ), offset = 1),
    #     column(3, wellPanel(
    #         uiOutput("fqUMIui"),
    #         shinyBS::bsTooltip(
    #             id = "fqUMIui", title = "List any unique molecular identifier (UMI) ranges to be extracted from the reads.",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         )
    #     )),
    #     column(3, wellPanel(
    #         uiOutput("fqCDNAui"),
    #         shinyBS::bsTooltip(
    #             id = "fqCDNAui", title = "List the cDNA range(s) in the reads to be mapped to the genome. May be one range (single-end) or two ranges (paired-end)",
    #             placement = "right", trigger = "hover", options = list(container = "body")
    #         )
    #     ))
    # ),
    # fluidRow(
    #     h4("Mapping/Reference options:", style = "padding-left: 20px;"),
    #     column(6, wellPanel(
    #         textInput(inputId = "STARpath", label = "Full path to the STAR index directory:", placeholder = "eg: /path/to/output"),
    #         shinyBS::bsTooltip(
    #             id = "STARpath", title = "The STAR index should be generated without splice junction database. Please remember to give the full path, as relative paths may not work.",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         ),
    #         textInput(inputId = "GTFpath", label = "Full path to the annotation GTF file:", placeholder = "eg: /path/to/output"),
    #         shinyBS::bsTooltip(
    #             id = "GTFpath", title = "Make sure the gene annotation matches the genome. Please remember to give the full path, as relative paths may not work.",
    #             placement = "bottom", trigger = "hover", options = list(container = "body")
    #         ),
    #         textInput(inputId = "STARparams", label = "Optional: Additional mapping parameters for STAR:", value = ""),
    #         shinyBS::bsTooltip(
    #             id = "STARparams", title = "You may list additional STAR mapping parameters. For instance, try trimming adapter sequenes using --clip3pAdapterSeq",
    #             placement = "top", trigger = "hover", options = list(container = "body")
    #         ),
    #         numericInput(inputId = "NUMadditionalFA", label = "Optional: Number of additional reference sequences:", value = 0, min = 0, step = 1),
    #         shinyBS::bsTooltip(
    #             id = "NUMadditionalFA", title = "Here you can give additional reference sequences PETRI-seq should map to but are not integrated in the STAR index. For instance, you could add the ERCC spike-in reference on the fly (eg. ERCC.fa).",
    #             placement = "top", trigger = "hover", options = list(container = "body")
    #         )
    #     )),
    #     column(6, wellPanel(
    #         p(strong("Additional references:")),
    #         uiOutput("refUI"),
    #         shinyBS::bsTooltip(
    #             id = "refUI", title = "Please remember to give the full paths, as relative paths may not work.",
    #             placement = "top", trigger = "hover", options = list(container = "body")
    #         )
    #     ))
)

load_yaml_panel <- tabPanel(
    "Load YAML",
    p(strong("You can load an existing YAML for PETRI-seq here")),
    textInput(inputId = "inYAML", label = "Full path to YAML file", placeholder = "eg: /path/to.yaml"),
    fileInput("file1", "...or choose YAML file to upload",
        accept = c(
            "text/tab-separated-values",
            "text/plain",
            ".yml",
            ".yaml"
        )
    ),
    actionButton(inputId = "LoadYAML", label = "Load YAML"),
    p(strong("Please click the Load button twice for all options to be properly set in shiny! Sorry for the inconvenience."))
)

save_yaml_panel <- tabPanel(
    "Generate YAML",
    p(strong("Your YAML is nearly ready...")),
    uiOutput("saveUI"),
    actionButton(inputId = "SaveYAML", label = "Save YAML!"),
    p(),
    p(strong("...or download your YAML file by clicking here:")),
    downloadButton("downloadData", "Download YAML")
)
