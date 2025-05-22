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

shortRead_to_Biostrings <- function(x) {
    QualityScaledDNAStringSet(
        sread(x),
        quality = quality(x)
    )
}

Biostrings_to_shortRead <- function(x) {
    ShortReadQ(
        sread = DNAStringSet(x),
        quality = quality(x)
    )
}
