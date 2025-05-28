suppressMessages({
    library(posDemux)
    library(glue)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]
bc_cutoff  <- as.integer(args[[2L]])

input_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")
freq_table <- read.table(file = input_frequency_table, header = TRUE, sep = "\t")

n_reads <- freq_table$frequency  |> sum()

plot_title <- glue("{sample}, n_reads={n_reads}") |> ggtitle()

freq_plot <- frequency_plot(freq_table, cutoff = bc_cutoff |> bc_to_frequency_cutoff(frequency_table = freq_table), type = "density", log_scale = TRUE) + plot_title
knee_plot <- knee_plot(freq_table, cutoff = bc_cutoff) + plot_title

ggsave(filename = "results/{sample}/{sample}_ReadsPerBC.pdf" |> glue(), plot = freq_plot)
ggsave(filename = "results/{sample}/{sample}_kneePlot.pdf" |> glue(), plot = knee_plot)
