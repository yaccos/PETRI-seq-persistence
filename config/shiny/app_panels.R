general_option_panel <- tabPanel(
    "General options",
    fluidRow(
        column(6, wellPanel(
            h4("Data and file input"),
            textInput("prefix", "Directory of input data", value = "resources/"),
            textInput("reference_genome", "Reference genome for alignment"),
            textOutput("path_reference_genome") |> h6(),
            textInput("reference_annotation", "Annotation GTF file"),
            textOutput("path_reference_annotation") |> h6(),
            textInput("suffix", "Common suffix for all read files", value = ""),
            textInput("forward_suffix", "Suffix for forward read files", value = ""),
            textInput("reverse_suffix", "Suffix for reverse read files", value = "")
        )),
        column(6, wellPanel(
            h4("Parameter options"),
            numericInput("chunk_size", "Chunk size for data streaming", value = 2e5, step = 10000),
            numericInput("bc_cutoff", "Number of barcode combinations to use", value = 1e4, step = 1000),
            shinyBS::bsTooltip("bc_cutoff",
             title = "Must be specified if you attempting to run the entire pipeline to completion. If you have not decided in advance, run the rule `determine_bc_cutoff` and follow the subsequent instructions.",
                placement = "bottom",
                trigger = "hover",
                options = list(container = "body")
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
    textInput("savePath", "Save YAML to", value = "config/config.yaml"),
    actionButton(inputId = "SaveYAML", label = "Save YAML"),
    p(),
    p(strong("...or download your YAML file by clicking here:")),
    downloadButton("downloadData", "Download YAML")
)
