suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
    library(ShortRead)
})

source("scripts/trim_helpers.R")

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]
threads  <- args[[2L]]

ShortRead:::.set_omp_threads(threads)

# sample  <- "random20000"

paired_input_file <- glue("results/{sample}/{sample}_QF_merged_R2_all_lanes.fastq")
input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")
input_bc_frame  <- glue("results/{sample}/{sample}_bc_frame.rds")
input_reads_to_keep  <- glue("results/{sample}/{sample}_selected_reads.txt")



message("Reading input sequences")
reverse_sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = paired_input_file, quality.scoring = "phred")  |> 
trim_sequence_names()

barcode_table  <- read.table(file = input_table_file, header = TRUE, row.names = NULL, sep = "\t")
bc_frame  <- readRDS(input_bc_frame)

reads_to_keep <- data.table::fread(input_reads_to_keep, nThread = 1L)[[1L]]

<<<<<<< HEAD:workflow/scripts/filter_and_trim_R2.R
bc1_stringset <- bc_frame$stringset$bc1
bc2_stringset <- bc_frame$stringset$bc2

message("Setting up FASTQ streams")

fq_chunk_size  <- as.integer(10^6)
fq_input_stream  <- FastqStreamer(paired_input_file, n = fq_chunk_size)
fq_output_stream  <- FastqOutput(output_file)


trim_sequence_names <- \(stringset) names(stringset) |>
    strsplit(" ") |>
    map_chr(1L)  |> 
    `names<-`(stringset, value=_)

process_chunk <- (chu)

handle_chunk <- function(chunk, outstream) {
    total_reads  <- length(chunk)
    chunk_ids  <- id(chunk)  |> sub(" .*$", "", _)
    keep_idx  <- chunk_ids %in% reads_to_keep  |> which()
    kept_reads  <- length(keep_idx)
    chunk_to_process <- chunk[keep_idx]
    chunk_to_process_ids  <- chunk_ids[keep_idx]
    chunk_barcode_match  <- match(chunk_to_process_ids, barcode_table$read)
    chunk_barcodes <- barcode_table[chunk_barcode_match,]
    chunk_bc1  <- chunk_barcodes$bc1
    chunk_bc1_assignment  <- split(seq_len(kept_reads), chunk_bc1)

    writeFastq(processed_reads, outstream)
=======
select_reads_from_cutoff <- function(filtered_res, bc_cutoff) {
    selected_freq_table <- freq_table[seq_len(cut_cutoff)]
    barcodes_to_keep <- selected_freq_table[bc_frame$bc_name]
    # barcode_table$rowID <- seq_len(nrow(barcode_table))
    common_rows <- dplyr::inner_join(barcode_table, barcodes_to_keep, by = names(barcodes_to_keep))
    kept_reads <- common_rows$read
    list(retained_reads = kept_reads, frequency_table = selected_freq_table)
>>>>>>> profile:workflow/scripts/filter_by_frequency.R
}






while (length(chunk  <- yield(fq_input_stream) > 0L)){

}

close(fq_input_stream)
close(fq_output_stream)

bc1_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc1", drop = TRUE]

message("Trimming BC1 and adapters from R2")




bc1_trimmed_R2_list <- imap(reverseComplement(bc1_stringset), trim_bc1_from_stringset)

bc1_trimmed_R2 <- rlang::exec(c, !!!unname(bc1_trimmed_R2_list) |> map("sequences"))
bc1_trim_count <- bc1_trimmed_R2_list |>
    map_int("trim_count") |>
    sum()
# Remove reads shorter than 16 nt
bc1_min_R2_length <- 16L
bc1_seq_too_short <- width(bc1_trimmed_R2) < bc1_min_R2_length
bc1_trimmed_R2 <- bc1_trimmed_R2[!bc1_seq_too_short]
trim_percentage <- bc1_trim_count / length(reverse_sequences_to_keep) * 100
n_sequences_too_short <- sum(bc1_seq_too_short)
percent_sequences_too_short <- n_sequences_too_short / length(reverse_sequences_to_keep) * 100
cat("Trimmed BC1 and its adjecent adapter from {bc1_trim_count} ({trim_percentage |> round(2L)}%) sequences\n\n" |> glue())
cat("After this trim, {n_sequences_too_short} ({percent_sequences_too_short |> round(2L)}%) \\
 sequences were removed because they became too short\n\n" |> glue())

message("Removing hairpins from R2")

trim_hairpins_from_stringset <- function(bc2_string, bc2_name) {
    # These sequences are *not* reverse compliments, but are aimed at hairpins
    bc3_to_bc2_adapter <- "GGTCCTTGGCTTCGC"
    bc2_to_bc1_adapter <- "CCTCCTACGCCAGA"
    combined_string <- xscat(bc3_to_bc2_adapter, bc2_string, bc2_to_bc1_adapter)

    sequences_with_bc <- bc2_selected_assigned_barcodes[names(bc1_trimmed_R2)] == bc2_name
    this_sequences <- bc1_trimmed_R2[sequences_with_bc]

    barcode_adapter_match <- Biostrings::vcountPattern(
        pattern = combined_string,
        subject = this_sequences,
        max.mismatch = 1L,
        with.indels = TRUE
    )

    hairpin_detected <- barcode_adapter_match != 0L
    this_trimmed_sequences <- this_sequences[!hairpin_detected]
    list(sequences = this_trimmed_sequences, trim_count = sum(hairpin_detected))
}

bc2_stringset <- bc_frame$stringset[[2]]
bc2_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc2", drop = TRUE]
names(bc2_selected_assigned_barcodes)  <- selected_assigned_barcodes$read
hairpin_trim <- imap(bc2_stringset, trim_hairpins_from_stringset)
hairpin_trimmed_R2 <- rlang::exec(c, !!!unname(hairpin_trim) |> map("sequences"))
hairpin_trim_count <- hairpin_trim |>
    map_int("trim_count") |>
    sum()
trim_percentage <- hairpin_trim_count / length(bc1_trimmed_R2) * 100
cat("Removed reads with hairpins from {hairpin_trim_count} ({trim_percentage |> round(2L)}%) sequences\n\n" |> glue())
writeQualityScaledXStringSet(hairpin_trimmed_R2, filepath = "results/{sample}/{sample}_2trim.fastq" |> glue(), compress = FALSE)
