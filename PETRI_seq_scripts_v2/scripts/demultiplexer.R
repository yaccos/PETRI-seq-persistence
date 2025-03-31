library(posDemux)
library(Biostrings)
library(purrr)
library(glue)
library(tibble)
library(ggplot2)
library(posDemux)

# These lines do not really belong here, but are used to debug and test the script in isolation
sample <- "random20000"
script_dir <- "/home/japet/Dokumenter/drug_presister_project/PETRI-seq-persistence/PETRI_seq_scripts_v2/demo/../scripts"

source(glue("{script_dir}/demultiplexer_helpers.R"))

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L

bc_frame <- tibble(bc_name = glue("bc{1:3}"))

barcode_file_prefix <- glue("{script_dir}/sc_barcodes_v2/")

barcode_file_name_main <- c("BC1_5p_anchor_v2.fa", "BC2_anchored.fa", "BC3_anchored.fa")


input_file <- glue("{sample}/{sample}_QF_UMI_L001_R1_001.fastq")

paired_input_file <- glue("{sample}/{sample}_QF_UMI_L001_R2_001.fastq")

output_table_file <- glue("{sample}_barcode_table.txt")

sequence_annotation <- c("B","A","B","A","B","A")

segment_lengths <- c(7L, 15L, 7L, 14L, 7L, NA_integer_)

min_sequence_length <- sum(head(segment_lengths, -1L))


bc_frame$filename <- paste0(barcode_file_prefix, barcode_file_name_main)

bc_frame$stringset <- map(bc_frame$filename, function(filepath) {
    message(glue("Reading filepath {filepath}"))
    raw_stringset <- Biostrings::readDNAStringSet(filepath = filepath)
    message(glue("Trimming away adapters in barcode file"))
    # The FASTA files contain the barcodes in addition to the adapters, so we must filter them out
    Biostrings::subseq(raw_stringset,start=1L, width = BARCODE_WIDTH)
    }
    )

names(bc_frame$stringset) <- bc_frame$bc_name


message("Reading input sequences")
sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = input_file, quality.scoring = "phred")
sequence_long_enough <- width(sequences) >= min_sequence_length
sequences <- sequences[sequence_long_enough]
paired_sequences <- Biostrings::readQualityScaledDNAStringSet(filepath = paired_input_file, quality.scoring = "phred")
paired_sequences <- paired_sequences[sequence_long_enough]

message("Starting demultiplexing")
demultiplex_res <- posDemux::combinatorial_demultiplex(sequences = sequences, barcodes = bc_frame$stringset  |> rev() , segments = sequence_annotation,
 segment_lengths = segment_lengths)

filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = ALLOWED_MISMATCHES)

res_table <- as.data.frame(filtered_res$demultiplex_res$assigned_barcodes)
res_table$read <- rownames(res_table)
res_table <- res_table[,c(4L,1L,2L,3L)]

write.table(res_table,
        output_table_file,
        sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = TRUE
    )

print(filtered_res$summary_res)

freq_table <- create_frequency_table(filtered_res$demultiplex_res$assigned_barcode)

freq_plot <- frequency_plot(freq_table)
knee_plot <- knee_plot(freq_table)

ggsave(filename = "{sample}_ReadsPerBC.pdf"  |> glue(), plot=freq_plot)
ggsave(filename= "{sample}_kneePlot.pdf" |> glue(), plot=knee_plot)
