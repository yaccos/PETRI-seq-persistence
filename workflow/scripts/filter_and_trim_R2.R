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

# HACK: We want to be able to control the number of threads used by ShortRead when parsing and writing FASTQ files
# According to the ShortRead documentation, this is the way to proceed even though it requires us to use an unexported function

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

barcodes_for_chunk(chunk, barcode){
    chunk_barcode_match  <- match(names(chunk), barcode_table$read)
    chunk_barcodes <- barcode_table[chunk_barcode_match,]
    chunk_barcodes[[barcode]] |> split(seq_along(chunk), _)
}


trim_bc1 <- function(chunk) {
    chunk_bc1 <- barcodes_for_chunk(chunk, "bc1")
    bc1_trimmed_R2_list  <- map2(chunk_bc1, bc1_stringset[names(chunk_bc1)] |> reverseComplement(),
     function(seq_idxs, stringset) trim_bc1_from_stringset(chunk_to_process[seq_idxs], stringset) )
    trim_count <- bc1_trimmed_R2_list |>
    map_int("trim_count") |>
    sum()
    bc1_trimmed_R2 <- rlang::exec(c, !!!unname(bc1_trimmed_R2_list) |> map("sequences"))

    # Remove reads shorter than 16 nt
    bc1_min_R2_length <- 16L
    bc1_seq_too_short <- width(bc1_trimmed_R2) < bc1_min_R2_length
    bc1_trimmed_R2 <- bc1_trimmed_R2[!bc1_seq_too_short]
    n_sequences_too_short <- sum(bc1_seq_too_short)

    list(chunk=bc1_trimmed_R2, statistics=list(bc1_trim_count = trim_count, n_bc1_too_short = n_sequences_too_short))
}

trim_hairpins <- function(chunk) {
    chunk_barcode_match  <- match(names(chunk), barcode_table$read)
    chunk_barcodes <- barcode_table[chunk_barcode_match,]
    chunk_barcodes <- 
    chunk_bc2_assignment <- split()

    bc2_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc2", drop = TRUE]
    names(bc2_selected_assigned_barcodes)  <- selected_assigned_barcodes$read
    hairpin_trim <- imap(bc2_stringset, trim_hairpins_from_stringset)
    hairpin_trimmed_R2 <- rlang::exec(c, !!!unname(hairpin_trim) |> map("sequences"))
    hairpin_trim_count <- hairpin_trim |>
    map_int("trim_count") |>
    sum()
}

filter_frequency(chunk) {
    total_reads  <- length(chunk)
    keep_idx  <- chunk_ids %in% reads_to_keep  |> which()
    n_kept_reads  <- length(keep_idx)
    list(chunk[keep_idx], statistics(total_reads=total_reads, kept_reads = n_kept_reads))
}

get_chunk_results(chunk){
    chunk_barcode_match  <- match(names(chunk), barcode_table$read)
    chunk_barcodes <- barcode_table[chunk_barcode_match,]
    chunk_bc1  <- chunk_barcodes$bc1
    chunk_bc1_assignment  <- split(seq_len(kept_reads), chunk_bc1)
    bc1_trimmed_R2_list  <- map2(chunk_bc1_assignment, bc1_stringset[names(chunk_bc1_assignment)] |> reverseComplement(),
     function(seq_idxs, stringset) trim_bc1_from_stringset(chunk_to_process[seq_idxs], stringset) )
    bc1_trim_count <- bc1_trimmed_R2_list |>
    map_int("trim_count") |>
    sum()
    bc1_trimmed_R2 <- rlang::exec(c, !!!unname(bc1_trimmed_R2_list) |> map("sequences"))

    # Remove reads shorter than 16 nt
    bc1_min_R2_length <- 16L
    bc1_seq_too_short <- width(bc1_trimmed_R2) < bc1_min_R2_length
    bc1_trimmed_R2 <- bc1_trimmed_R2[!bc1_seq_too_short]
    n_sequences_too_short <- sum(bc1_seq_too_short)

    chunk_bc2_assignment <- split()

    bc2_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc2", drop = TRUE]
    names(bc2_selected_assigned_barcodes)  <- selected_assigned_barcodes$read
    hairpin_trim <- imap(bc2_stringset, trim_hairpins_from_stringset)
    hairpin_trimmed_R2 <- rlang::exec(c, !!!unname(hairpin_trim) |> map("sequences"))
    hairpin_trim_count <- hairpin_trim |>
    map_int("trim_count") |>
    sum()
}

handle_chunk <- function(chunk, outstream) {
    chunk_ids  <- id(chunk)  |> sub(" .*$", "", _)
    chunk |>
     ShortRead_to_Biostrings()  |>
     `names<-(chunk_ids)`() |>  
     filter_chain(filter_frequency, trim_bc1, trim_hairpins)
    

    writeFastq(processed_reads, outstream)
}


while (length(chunk  <- yield(fq_input_stream) > 0L)){
    
}

close(fq_input_stream)
close(fq_output_stream)

bc1_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc1", drop = TRUE]

message("Trimming BC1 and adapters from R2")



trim_percentage <- bc1_trim_count / length(reverse_sequences_to_keep) * 100

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
