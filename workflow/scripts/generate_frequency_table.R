suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(tibble)
    library(data.table)
})

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

input_table  <- "results/{sample}/{sample}_barcode_table.txt"

output_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

barcode_table <- fread(input_table, sep = "\t", header = TRUE)

assigned_barcode  <- matrix(barcode_table[,c("bc1","bc2","bc3")], dimnames = list(barcode_table$read, names(barcode_table)))

freq_table <- create_frequency_table(assigned_barcode)

write.table(x = freq_table, file = output_frequency_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
