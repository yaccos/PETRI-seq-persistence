library(posDemux)
library(Biostrings)
library(purrr)
library(glue)
library(tibble)
library(dplyr)
library(ggplot2)
library(Rsamtools)

# These lines do not really belong here, but are used to debug and test the script in isolation
sample <- "random20000"
script_dir <- "/home/japet/Dokumenter/drug_presister_project/PETRI-seq-persistence/PETRI_seq_scripts_v2/demo/../scripts"

sam_file  <- glue("{sample}/{sample}_bwa_sam.sam")

