library(Biostrings)
library(glue)
library(ggplot2)
library(dplyr)
library(purrr)


# args <- commandArgs(trailingOnly = TRUE)

# sample  <- args[1L]
# lane  <- args[2L]

# Only for debugging, remove before production
sample <- "random20000"
lane <- "1"

min_overlap_to_trim <- 8L
min_length_after_trim <- 75L
forward_length_to_keep <- 58L

forward_fastq <- glue("results/{sample}/{sample}_QF_L00{lane}_R1_001.fastq")
reverse_fastq <- glue("results/{sample}/{sample}_QF_L00{lane}_R2_001.fastq")

trim_sequence_names <- \(stringset) names(stringset) |>
    strsplit(" ") |>
    map_chr(1L)  |> 
    `names<-`(stringset, value=_)

forward_reads <- Biostrings::readQualityScaledDNAStringSet(filepath = forward_fastq, quality.scoring = "phred")  |>
 trim_sequence_names()
reverse_reads <- Biostrings::readQualityScaledDNAStringSet(filepath = reverse_fastq, quality.scoring = "phred")  |>
 trim_sequence_names()

overlap_alignment <- Biostrings::pairwiseAlignment(
    forward_reads,
    reverse_reads |> reverseComplement(),
    type = "overlap"
)

n_overlap_bases <- nchar(overlap_alignment)

overlap_frame <- data.frame(
    overlap = n_overlap_bases,
    forward_length = lengths(forward_reads),
    reverse_length = lengths(reverse_reads),
    seq_name = names(forward_reads)
) |>
    mutate(combined_length = forward_length + reverse_length, to_trim = overlap >= min_overlap_to_trim) |>
    # We have to subtract the overlap twice because it occurs on both reads
    mutate(combined_length_after_trim = ifelse(to_trim, combined_length - 2L * overlap, combined_length)) |>
    mutate(read_too_short = combined_length_after_trim < min_length_after_trim)  |> 
    mutate(reverse_trim_pos = ifelse(to_trim, reverse_length, reverse_length - overlap))

frame_to_keep  <- filter(overlap_frame, !read_too_short)

kept_forward_fastq  <- forward_reads  |> _[frame_to_keep$seq_name]  |> subseq(start = 1L, width = forward_length_to_keep)

cat("Results of overlap trimming:\n")
n_reads <- nrow(overlap_frame)
glue("Total number of reads: {n_reads}") |>
    cat("\n")
n_removed <- overlap_frame$read_too_short |> sum()
removed_percentage <- n_removed / n_reads * 100
"{n_removed} ({removed_percentage |> round(2L)}%) reads were removed because the paired reads had too much overlap" |>
    glue() |>
    cat("\n")
n_trimmed <- (overlap_frame$to_trim * !overlap_frame$read_too_short)  |> sum()
trimmed_percentage <- n_trimmed / n_reads * 100
"{n_trimmed} ({trimmed_percentage |> round(2L)}%) of the reads pairs were trimmed for the overlap" |>
    glue() |>
    cat("\n")


trim_histogram <- ggplot(overlap_frame, mapping = aes(x = overlap)) +
    stat_count() +
    xlab("Length of overlap") +
    ylab("Number of read pairs")

ggsave(filename = "results/{sample}_overlap_histogram.pdf" |> glue(), plot = trim_histogram)

