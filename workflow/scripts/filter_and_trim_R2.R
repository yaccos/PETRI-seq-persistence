suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
    library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)

# We must parse the arguments before sourcing the scripts because the scripts need to read 
# files of their own
sample  <- args[[1L]]
threads  <- args[[2L]]

source("scripts/trim_filters.R")


# HACK: We want to be able to control the number of threads used by ShortRead when parsing and writing FASTQ files
# According to the ShortRead documentation, this is the way to proceed even though it requires us to use an unexported function
threads <- .Call(ShortRead:::.set_omp_threads, 1L) 


paired_input_file <- glue("results/{sample}/{sample}_QF_merged_R2_all_lanes.fastq")
input_reads_to_keep  <- glue("results/{sample}/{sample}_selected_reads.txt")
output_file  <- glue("results/{sample}/{sample}_2trim.fastq")

reads_to_keep <- data.table::fread(input_reads_to_keep, nThread = 1L)[[1L]]


message("Setting up FASTQ streams")

fq_chunk_size  <- as.integer(10^6)
fq_input_stream  <- FastqStreamer(paired_input_file, n = fq_chunk_size)

processing_chain <- filter_chain(filter_frequency, trim_bc1, filter_bc1, trim_hairpins, kept_reads_count)

report_progress <- function(counts) {
    iteration_gap  <- as.integer(10^6)
    with(as.list(counts), {
        if(total_reads %% iteration_gap == 0L) {
            message(glue("Processed {total_reads} reads, kept {kept_reads} so far..."))    
        }
    }
    )
}

summarize_results <- function(counts) {
    with(as.list(counts), {
    message(glue("Completed processing {total_reads} reads"))
    message(glue("of which {selected_reads} reads ({selected_reads/total_reads*100}%) were selected from barcode frequency"))
    message(glue("Trimmed BC1 from {bc1_trim_count} reads"))
    message(glue("Removed {n_bc1_too_short} reads that were too short after trimming"))
    message(glue("Removed {hairpin_filter_count} reads with hairpins"))
    message(glue("A total of {kept_reads} ({kept_reads/selected_reads*100}%) of the selected reads were written to the output FASTQ file"))
    }
    )
}

# Initializes all counts to zero by calling the filtering chain on an empty input stringset
message("Initializing counts")
counts <- processing_chain(empty_chunk) |> _$counts
message("Starting streaming")
while ((chunk  <- yield(fq_input_stream))  |> length() > 0L) {
    chunk_ids  <- id(chunk)  |> sub(" .*$", "", x=_)
    chain_results <- chunk |>
    ShortRead_to_Biostrings()  |>
    `names<-`(chunk_ids) |>  
    processing_chain()
    processed_reads <- chain_results$chunk
    counts  <- counts + chain_results$counts
    writeQualityScaledXStringSet(processed_reads, output_file, append = TRUE)
    report_progress(counts)
}

close(fq_input_stream)

summarize_results(counts)
