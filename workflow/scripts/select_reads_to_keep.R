suppressMessages({
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
    library(data.table)
})

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

# sample  <- "random20000"

# bc_cutoff <- 7000L

bc_cutoff  <- as.integer(args[[2L]])

input_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

bc_names  <- glue("bc{1:3}")

log_progress("Reading barcode table")

barcode_table  <- fread(file = input_table_file, header = TRUE, sep = "\t")

log_progress("Reading frequency table")

freq_table <- read.table(file = input_frequency_table, header = TRUE, row.names = NULL, sep = "\t")

selected_freq_table <- freq_table[seq_len(bc_cutoff),]

log_progress("Selecting reads to keep")
common_rows <- posDemux::row_match(barcode_table, selected_freq_table)

pasted_celltag <- barcode_table[common_rows, do.call(paste, c(.SD, sep = "_")), .SDcols = bc_names]
selected_barcode_table  <- data.table(read=barcode_table$read[common_rows], UMI=barcode_table$UMI[common_rows],
 celltag=pasted_celltag)

kept_reads <- barcode_table$read[common_rows]

log_progress("Writing processed frequency table")
write.table(
    x = selected_freq_table, file = "results/{sample}/{sample}_selected_frequency_table.txt" |> glue(),
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

log_progress("Writing processed barcode table")
data.table::fwrite(x = selected_barcode_table, file = "results/{sample}/{sample}_selected_barcode_table.txt" |> glue(), sep="\t", nThread = 1L)

log_progress("DONE")

n_selected <- sum(common_rows)
n_read  <- nrow(barcode_table)
selection_percentage <- n_selected/n_read*100
glue("Selected {n_selected} of {n_read} reads ({round(selection_percentage, 2L)}%) based on barcode frequency")
