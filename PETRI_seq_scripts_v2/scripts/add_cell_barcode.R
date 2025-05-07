# This script is intended to add a tag CB containing the cell barcode and a tag BX containing the UMI
# to a BAM file processed through featureCounts
suppressMessages({
    library(GenomicAlignments)
    library(glue)
    library(rtracklayer)
})
args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]

bam_file <- glue("results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam")
barcode_index_file <- glue("results/{sample}_barcode_table.txt")

annotated_alignment <- import(bam_file, format = "BAM", use.names = TRUE, param = ScanBamParam(tag = c("X0", "XA", "NM", "MD", "XT")))
barcode_index <- read.table(barcode_index_file, header = TRUE)
barcode_index$celltags <- rlang::exec(paste, !!!barcode_index[c("bc1", "bc2", "bc3")], sep = "_")
lookup_matches <- match(names(annotated_alignment), barcode_index$read)
mcols(annotated_alignment)$CB <- barcode_index$celltags[lookup_matches]
mcols(annotated_alignment)$BX <- barcode_index$UMI[lookup_matches]

output_file <- glue("results/{sample}_FC/{sample}_sorted.bam.featureCounts_with_celltag.bam")
export(annotated_alignment, output_file)
