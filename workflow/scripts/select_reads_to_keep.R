suppressMessages({
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

# sample  <- "random20000"

# bc_cutoff <- 7000L

bc_cutoff  <- as.integer(args[[2L]])

input_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

bc_names  <- glue("bc{1:3}")

barcode_table  <- read.table(file = input_table_file, header = TRUE, row.names = NULL, sep = "\t")

trim_sequence_names <- \(stringset) names(stringset) |>
    strsplit(" ") |>
    map_chr(1L)  |> 
    `names<-`(stringset, value=_)

freq_table <- read.table(file = input_frequency_table, header = TRUE, row.names = NULL, sep = "\t")

selected_freq_table <- freq_table[seq_len(bc_cutoff)]
barcodes_to_keep <- selected_freq_table[bc_names]
common_rows <- dplyr::inner_join(barcode_table, barcodes_to_keep, by = names(barcodes_to_keep))
kept_reads <- common_rows$read

write.table(
    x = selected_freq_table, file = "results/{sample}/{sample}_selected_frequency_table.txt" |> glue(),
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

data.table::fwrite(x = kept_reads |> list(), file = "results/{sample}/{sample}_selected_reads.txt" |> glue(), nThread = 1L)
