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

barcode_names  <- c("bc1","bc2","bc3")

input_table  <- glue("results/{sample}/{sample}_barcode_table.txt")

output_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

barcode_table <- fread(input_table, sep = "\t", header = TRUE)

# Usually, this matrix contains the read identifiers as rownames, but it is not necessary here
# The column names are essential for creating the frequency table, but they are automatically inferred from
# the names of the barcodes
assigned_barcode  <- as.matrix(barcode_table[, ..barcode_names])

freq_table <- create_frequency_table(assigned_barcode)

write.table(x = freq_table, file = output_frequency_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
