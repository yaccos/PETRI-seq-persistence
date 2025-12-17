log_file = snakemake@log[[1]]
log_handle  <- file(log_file, open = "w")
sink(log_handle, append = TRUE, type = "output")
sink(log_handle, append = TRUE, type = "message")

suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(tibble)
    library(ShortRead)
})

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

bc_frame <- tibble(bc_name = glue("bc{1:3}"))
bc_frame$filename <- snakemake@input[bc_frame$bc_name] |> as.vector()

input_file <- snakemake@input[["fastq"]]

output_table_file <- snakemake@output[["barcode_table"]]
output_bc_frame  <- snakemake@output[["bc_frame"]]
output_freq_table  <- snakemake@output[["freq_table"]]
chunk_size <- snakemake@params[["chunk_size"]]

sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")

segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)

bc_frame$stringset <- map(bc_frame$filename, function(filepath) {
    glue("Reading filepath {filepath}")  |>  log_progress()
    raw_stringset <- Biostrings::readDNAStringSet(filepath = filepath)
    glue("Trimming away adapters in barcode file")  |> log_progress()
    # The FASTA files contain the barcodes in addition to the adapters, so we must filter them out
    Biostrings::subseq(raw_stringset, start = 1L, width = BARCODE_WIDTH)
})

names(bc_frame$stringset) <- bc_frame$bc_name

callbacks <- streaming_callbacks(input_file = input_file,
                                 output_table_file = output_table_file,
                                 chunk_size = chunk_size,
                                 verbose = TRUE)

streaming_res <- rlang::exec(streaming_demultiplex, !!! callbacks,
                             barcodes=bc_frame$stringset |> rev(), allowed_mismatches = 1L,
            segments = sequence_annotation, segment_lengths = segment_lengths)

freq_table  <- streaming_res$freq_table
log_progress("Writing frequency table...")
write.table(x = freq_table, file = output_freq_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
log_progress("DONE")
cat("\n")
print(streaming_res$summary_res)

saveRDS(object = bc_frame, file = output_bc_frame, compress = FALSE)

sink(type="message")
sink(type="output")
close(log_handle)
