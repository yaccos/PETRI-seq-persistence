suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(tibble)
    library(ShortRead)
    library(data.table)
})

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

chunk_size <- args[[2L]] |> as.integer()

bc_frame <- tibble(bc_name = glue("bc{1:3}"))

barcode_file_prefix <- glue("data/sc_barcodes_v2/")

barcode_file_name_main <- c("BC1_5p_anchor_v2.fa", "BC2_anchored.fa", "BC3_anchored.fa")

input_file <- glue("results/{sample}/{sample}_QF_R1_all_lanes.fastq")

output_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

output_bc_frame  <- glue("results/{sample}/{sample}_bc_frame.rds")

output_freq_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")

segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)


bc_frame$filename <- paste0(barcode_file_prefix, barcode_file_name_main)

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
