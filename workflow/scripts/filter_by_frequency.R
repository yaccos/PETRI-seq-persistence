# Prologue
suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

# sample  <- "random20000"

# bc_cutoff <- 7000L

bc_cutoff  <- as.integer(args[[2L]])

paired_input_file <- glue("results/{sample}/{sample}_QF_merged_R2_all_lanes.fastq")

input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

input_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")

input_bc_frame  <- glue("results/{sample}/{sample}_bc_frame.rds")

trim_sequence_names <- \(stringset) names(stringset) |>
    strsplit(" ") |>
    map_chr(1L)  |> 
    `names<-`(stringset, value=_)


message("Reading input sequences")
reverse_sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = paired_input_file, quality.scoring = "phred")  |> 
trim_sequence_names()

freq_table <- read.table(file = input_frequency_table, header = TRUE, row.names = NULL, sep = "\t")
barcode_table  <- read.table(file = input_table_file, header = TRUE, row.names = NULL, sep = "\t")
bc_frame  <- readRDS(input_bc_frame)

# Main context
# bc_cutoff <- posDemux::interactive_bc_cutoff(freq_table) |> print()

select_reads_from_cutoff <- function(filtered_res, bc_cutoff) {
    # The frequency table is already sorted, so we don't have to do anything fancy here
    selected_freq_table <- freq_table[seq_len(bc_cutoff)]
    barcodes_to_keep <- selected_freq_table[bc_frame$bc_name]
    # barcode_table$rowID <- seq_len(nrow(barcode_table))
    common_rows <- dplyr::inner_join(barcode_table, barcodes_to_keep, by = names(barcodes_to_keep))
    kept_reads <- common_rows$read
    list(retained_reads = kept_reads, frequency_table = selected_freq_table)
}

selection_res <- select_reads_from_cutoff(filtered_res, bc_cutoff)

reads_to_keep <- selection_res$retained_reads
selected_frequency_table <- selection_res$frequency_table

write.table(
    x = selected_frequency_table, file = "results/{sample}/{sample}_selected_frequency_table.txt" |> glue(),
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

reverse_sequences_to_keep <- reverse_sequences[reads_to_keep]
selected_assigned_barcodes <- barcode_table  |> filter(read %in% reads_to_keep)

bc1_stringset <- bc_frame$stringset[[1]]
bc1_selected_assigned_barcodes <- selected_assigned_barcodes[, "bc1", drop = TRUE]

message("Trimming BC1 and adapters from R2")

trim_bc1_from_stringset <- function(bc_string, bc_name) {
    # This is the reverse compliment of the adapter found between BC1 and BC2
    trimming_adapter_sequence <- "TCTGGCGTAGGAGG"
    this_R2_sequences <- reverse_sequences_to_keep[bc1_selected_assigned_barcodes == bc_name]

    adapter_match <- Biostrings::vmatchPattern(
        pattern = trimming_adapter_sequence,
        subject = this_R2_sequences,
        max.mismatch = 1L,
        with.indels = FALSE
    )

    barcode_adapter_match <- Biostrings::vmatchPattern(
        pattern = xscat(bc_string, trimming_adapter_sequence),
        subject = this_R2_sequences,
        max.mismatch = 2L,
        with.indels = FALSE
    )

    extract_first_match <- function(match_object) {
        match_lengths <- lengths(match_object)
        res <- rep(NA_integer_, length(this_R2_sequences))
        sequences_with_match <- match_lengths > 0L
        res[sequences_with_match] <- map_int(match_object[sequences_with_match] |> startIndex(), 1L)
        res
    }

    match_results <- map(list(adapter_match, barcode_adapter_match), extract_first_match)
    first_combined_match <- rlang::exec(pmin, !!!match_results, na.rm = TRUE)

    this_trimmed_sequences <- this_R2_sequences
    sequences_to_trim <- !is.na(first_combined_match)
    this_trimmed_sequences[sequences_to_trim] <- subseq(this_R2_sequences[sequences_to_trim], end = first_combined_match[sequences_to_trim] - 1L)
    list(sequences = this_trimmed_sequences, trim_count = sum(sequences_to_trim))
}


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
