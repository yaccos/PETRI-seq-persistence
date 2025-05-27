source("scripts/trim_helpers.R")

input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")
input_bc_frame  <- glue("results/{sample}/{sample}_bc_frame.rds")

barcode_table  <- read.table(file = input_table_file, header = TRUE, row.names = NULL, sep = "\t")
bc_frame  <- readRDS(input_bc_frame)


bc1_stringset <- bc_frame$stringset$bc1
bc2_stringset <- bc_frame$stringset$bc2

barcodes_for_chunk <- function(chunk, barcode) {
    chunk_barcode_match  <- match(names(chunk), barcode_table$read)
    chunk_barcodes <- barcode_table[chunk_barcode_match,]
    chunk_barcodes[[barcode]] |> split(seq_along(chunk), f=_)
}

empty_chunk <- character() |> DNAStringSet()

concatenate_reads <- function(read_list){
    # There is a tricky edge case when the input chunk is empty
    if (length(read_list) == 0L) {
        empty_chunk
    } else {
        rlang::exec(c, !!!unname(read_list) |> map("sequences"))
    }
}

trim_bc1 <- function(chunk) {
    chunk_bc1 <- barcodes_for_chunk(chunk, "bc1")
    bc1_trimmed_R2_list  <- map2(chunk_bc1, bc1_stringset[names(chunk_bc1)] |> reverseComplement(),
     function(seq_idxs, stringset) trim_bc1_from_stringset(chunk[seq_idxs], stringset) )
    trim_count <- bc1_trimmed_R2_list |>
    map_int("trim_count") |>
    sum()
    trimmed_chunk <- concatenate_reads(bc1_trimmed_R2_list)
    
    list(chunk = trimmed_chunk, counts = c(bc1_trim_count = trim_count))
}

filter_bc1 <- function(chunk) {
    # Remove reads shorter than 16 nt
    bc1_min_R2_length <- 16L
    bc1_seq_too_short <- width(chunk) < bc1_min_R2_length
    filtered_chunk <- chunk[!bc1_seq_too_short]
    n_sequences_too_short <- sum(bc1_seq_too_short)

    list(chunk = filtered_chunk, counts = c(n_bc1_too_short = n_sequences_too_short))
}


trim_hairpins <- function(chunk) {
    chunk_bc2 <- barcodes_for_chunk(chunk, "bc2")

    hairpin_filter <-  map2(chunk_bc2, bc2_stringset[names(chunk_bc2)] |> reverseComplement(),
     function(seq_idxs, stringset) filter_hairpins_from_stringset(chunk[seq_idxs], stringset) )

    filtered_chunk <- concatenate_reads(hairpin_filter)
    hairpin_filter_count <- hairpin_filter |>
    map_int("trim_count") |>
    sum()
    list(chunk=filtered_chunk, counts = c(hairpin_filter_count = hairpin_filter_count))
}

filter_frequency <- function(chunk) {
    total_reads  <- length(chunk)
    keep_idx  <- names(chunk) %in% reads_to_keep  |> which()
    n_kept_reads  <- length(keep_idx)
    list(chunk=chunk[keep_idx], counts = c(total_reads=total_reads, selected_reads = n_kept_reads))
}

kept_reads_count <- function(chunk) {
    kept_reads = length(chunk)
    list(chunk=chunk, counts = c(kept_reads=kept_reads))
}
