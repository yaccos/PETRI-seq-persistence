trim_bc1_from_stringset <- function(chunk, bc_string) {
    # chunk is assumed to be a QualityScaledDNAStringSet
    # This is the reverse compliment of the adapter found between BC1 and BC2
    trimming_adapter_sequence <- "TCTGGCGTAGGAGG"

    adapter_match <- Biostrings::vmatchPattern(
        pattern = trimming_adapter_sequence,
        subject = chunk,
        max.mismatch = 1L,
        with.indels = FALSE
    )

    barcode_adapter_match <- Biostrings::vmatchPattern(
        pattern = xscat(bc_string, trimming_adapter_sequence),
        subject = chunk,
        max.mismatch = 2L,
        with.indels = FALSE
    )

    extract_first_match <- function(match_object) {
        match_lengths <- lengths(match_object)
        res <- rep(NA_integer_, length(chunk))
        sequences_with_match <- match_lengths > 0L
        res[sequences_with_match] <- map_int(match_object[sequences_with_match] |> startIndex(), 1L)
        res
    }

    match_results <- map(list(adapter_match, barcode_adapter_match), extract_first_match)
    first_combined_match <- rlang::exec(pmin, !!!match_results, na.rm = TRUE)

    trimmed_sequences <- chunk
    sequences_to_trim <- !is.na(first_combined_match)
    trimmed_sequences[sequences_to_trim] <- subseq(chunk[sequences_to_trim], end = first_combined_match[sequences_to_trim] - 1L)
    list(sequences = trimmed_sequences, trim_count = sum(sequences_to_trim))
}

trim_hairpins_from_stringset <- function(chunk, bc2_string, bc2_name) {
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

filter_chain <- function(chunk, ...) {
    # This is a special function which accepts a chunk and any given number of functions of the form f(chunk) -> list(chunk, statistics)
    # and chains the operations the functions together
    statistics  <- list()
    for (f in ...){
        f_res  <- f(chunk)
        statistics  <- c(statistics, f_res$statistics)
        chunk <- f_res$chunk
    }
    list(chunk =  chunk, statistics = statistics)
}

ShortRead_to_Biostrings <- function(x) {
    QualityScaledDNAStringSet(
        sread(x),
        quality = quality(x)
    )
}

Biostrings_to_ShortRead <- function(x) {
    ShortReadQ(
        sread = DNAStringSet(x),
        quality = quality(x)
    )
}
