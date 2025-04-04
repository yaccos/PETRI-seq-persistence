library(posDemux)
library(Biostrings)
library(purrr)
library(glue)
library(tibble)
library(dplyr)
library(ggplot2)
library(posDemux)

# These lines do not really belong here, but are used to debug and test the script in isolation
sample <- "random20000"
script_dir <- "/home/japet/Dokumenter/drug_presister_project/PETRI-seq-persistence/PETRI_seq_scripts_v2/demo/../scripts"

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

bc_frame <- tibble(bc_name = glue("bc{1:3}"))

barcode_file_prefix <- glue("{script_dir}/sc_barcodes_v2/")

barcode_file_name_main <- c("BC1_5p_anchor_v2.fa", "BC2_anchored.fa", "BC3_anchored.fa")


input_file <- glue("{sample}/{sample}_QF_UMI_L001_R1_001.fastq")

paired_input_file <- glue("{sample}/{sample}_QF_UMI_L001_R2_001.fastq")

output_table_file <- glue("{sample}_barcode_table.txt")

sequence_annotation <- c("B", "A", "B", "A", "B", "A")

segment_lengths <- c(7L, 15L, 7L, 14L, 7L, NA_integer_)

min_sequence_length <- sum(head(segment_lengths, -1L))


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
forward_sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = input_file, quality.scoring = "phred")
sequence_long_enough <- width(forward_sequences) >= min_sequence_length
forward_sequences <- forward_sequences[sequence_long_enough]
reverse_sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = paired_input_file, quality.scoring = "phred")
reverse_sequences <- reverse_sequences[sequence_long_enough]

message("Starting demultiplexing")
demultiplex_res <- posDemux::combinatorial_demultiplex(
    sequences = forward_sequences, barcodes = bc_frame$stringset |> rev(), segments = sequence_annotation,
    segment_lengths = segment_lengths
)

filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = ALLOWED_MISMATCHES)

res_table <- as.data.frame(filtered_res$demultiplex_res$assigned_barcodes)
res_table$read <- rownames(res_table)
res_table <- res_table[, c(4L, 1L, 2L, 3L)]

write.table(res_table,
    output_table_file,
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

print(filtered_res$summary_res)

freq_table <- create_frequency_table(filtered_res$demultiplex_res$assigned_barcode)

write.table(x = freq_table, file = "{sample}_frequency_table.txt" |> glue(), quote = FALSE, sep = "\t", row.names = FALSE)

# bc_cutoff <- posDemux::interactive_bc_cutoff(freq_table) |> print()
bc_cutoff <- 7000L


freq_plot <- frequency_plot(freq_table, cutoff = bc_cutoff |> bc_to_frequency_cutoff(frequency_table = freq_table))
knee_plot <- knee_plot(freq_table, cutoff = bc_cutoff)

ggsave(filename = "{sample}_ReadsPerBC.pdf" |> glue(), plot = freq_plot)
ggsave(filename = "{sample}_kneePlot.pdf" |> glue(), plot = knee_plot)


select_reads_from_cutoff <- function(filtered_res, bc_cutoff) {
    freq_table <- create_frequency_table(filtered_res$demultiplex_res$assigned_barcodes)
    barcode_table <- filtered_res$demultiplex_res$assigned_barcodes |> as.data.frame()
    selected_freq_table <- freq_table[freq_table$cumulative_frequency <= bc_cutoff, ]
    barcodes_to_keep <- selected_freq_table[colnames(barcode_table)]
    barcode_table$rowID <- seq_len(nrow(barcode_table))
    common_rows <- dplyr::inner_join(barcode_table, barcodes_to_keep, by = names(barcodes_to_keep))
    kept_reads <- common_rows$rowID
    retained_reads <- rep_along(filtered_res$retained, FALSE)
    retained_reads[kept_reads] <- TRUE
    # The resulting reads must have a barcode frequent enough AND have a number of mismatches below the
    # selected threshold
    list(retained_reads = retained_reads & filtered_res$retained, frequency_table = selected_freq_table)
}

selection_res <- select_reads_from_cutoff(filtered_res, bc_cutoff)

reads_to_keep <- selection_res$retained_reads
selected_frequency_table <- selection_res$frequency_table

write.table(
    x = selected_frequency_table, file = "{sample}_selected_frequency_table.txt" |> glue(),
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

forward_sequences_to_keep  <- forward_sequences[reads_to_keep]
reverse_sequences_to_keep <- reverse_sequences[reads_to_keep]
selected_assigned_barcodes  <- demultiplex_res$assigned_barcodes[reads_to_keep, ]

bc1_stringset <- bc_frame$stringset[[1]]
bc1_selected_assigned_barcodes  <- selected_assigned_barcodes[, "bc1", drop=TRUE]

trim_bc1_from_stringset  <- function(bc_string, bc_name) {
    # This is the reverse compliment of the adapter found between BC1 and BC2
    trimming_adapter_sequence <- "TCTGGCGTAGGAGG"
    this_sequences <- reverse_sequences_to_keep[bc1_selected_assigned_barcodes == bc_name]

    adapter_match  <- Biostrings::vmatchPattern(pattern = trimming_adapter_sequence,
     subject = this_sequences,
     max.mismatch = 1L,
     with.indels = FALSE
    )

    barcode_adapter_match  <- Biostrings::vmatchPattern(
        pattern = xscat(bc_string, trimming_adapter_sequence),
     subject = this_sequences,
     max.mismatch = 2L,
     with.indels = FALSE
    )

    extract_first_match  <- function(match_object) {
        match_lengths <- lengths(match_object)
        res  <- rep(NA_integer_, length(this_sequences))
        sequences_with_match  <- match_lengths > 0L
        res[sequences_with_match] <- map_int(match_object[sequences_with_match] |> startIndex(), 1L)
        res
    }

    match_results <- map(list(adapter_match, barcode_adapter_match), extract_first_match)
    first_combined_match  <- rlang::exec(pmin, !!! match_results, na.rm=TRUE)

    this_trimmed_sequences  <- this_sequences
    sequences_to_trim  <- !is.na(first_combined_match)
    this_trimmed_sequences[sequences_to_trim] <- subseq(this_sequences[sequences_to_trim], end = first_combined_match[sequences_to_trim] - 1L)
    this_trimmed_sequences
}


bc1_trimmed_R2_list <- imap(reverseComplement(bc1_stringset), trim_bc1_from_stringset)

bc1_trimmed_R2  <- rlang::exec(c, !!! unname(bc1_trimmed_R2_list))
# Remove reads shorter than 16 nt
bc1_min_R2_length  <- 16L
bc1_seq_too_short  <- width(bc1_trimmed_R2) < bc1_min_R2_length
bc1_trimmed_R1 <- forward_sequences_to_keep[!bc1_seq_too_short]
bc1_trimmed_R2 <- bc1_trimmed_R2[!bc1_seq_too_short]
writeQualityScaledXStringSet(bc1_trimmed_R2, filepath = "{sample}/{sample}_R2_trimmed.fastq" |> glue(), compress = FALSE)
