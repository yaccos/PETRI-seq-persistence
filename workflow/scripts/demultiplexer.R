suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(tibble)
})

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

bc_frame <- tibble(bc_name = glue("bc{1:3}"))

barcode_file_prefix <- glue("data/sc_barcodes_v2/")

barcode_file_name_main <- c("BC1_5p_anchor_v2.fa", "BC2_anchored.fa", "BC3_anchored.fa")

input_file <- glue("results/{sample}/{sample}_QF_merged_R1_all_lanes.fastq")

output_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

output_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

output_bc_frame  <- glue("results/{sample}/{sample}_bc_frame.rds")

sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")

segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)

min_sequence_length <- sum(head(segment_lengths, -1L))

trim_sequence_names <- \(stringset) names(stringset) |>
    strsplit(" ") |>
    map_chr(1L)  |> 
    `names<-`(stringset, value=_)


bc_frame$filename <- paste0(barcode_file_prefix, barcode_file_name_main)

bc_frame$stringset <- map(bc_frame$filename, function(filepath) {
    message(glue("Reading filepath {filepath}"))
    raw_stringset <- Biostrings::readDNAStringSet(filepath = filepath)
    message(glue("Trimming away adapters in barcode file"))
    # The FASTA files contain the barcodes in addition to the adapters, so we must filter them out
    Biostrings::subseq(raw_stringset, start = 1L, width = BARCODE_WIDTH)
})

names(bc_frame$stringset) <- bc_frame$bc_name

message("Reading input sequences")
forward_sequences <- Biostrings::readDNAStringSet(filepath = input_file, format = "fastq")  |> 
trim_sequence_names()

sequence_long_enough <- width(forward_sequences) >= min_sequence_length
forward_sequences <- forward_sequences[sequence_long_enough]

message("Starting demultiplexing")
demultiplex_res <- posDemux::combinatorial_demultiplex(
    sequences = forward_sequences, barcodes = bc_frame$stringset |> rev(), segments = sequence_annotation,
    segment_lengths = segment_lengths
)

filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = ALLOWED_MISMATCHES)

res_table <- as.data.frame(filtered_res$demultiplex_res$assigned_barcodes)
res_table$UMI <- filtered_res$demultiplex_res$payload$UMI |> as.character()
res_table$read <- rownames(res_table)
res_table <- res_table[, c("read", "UMI", "bc3", "bc2", "bc1")]

# Export of results

write.table(res_table,
    output_table_file,
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

print(filtered_res$summary_res)

freq_table <- create_frequency_table(filtered_res$demultiplex_res$assigned_barcode)

write.table(x = freq_table, file = output_frequency_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

saveRDS(object = bc_frame, file = output_bc_frame, compress = FALSE)
