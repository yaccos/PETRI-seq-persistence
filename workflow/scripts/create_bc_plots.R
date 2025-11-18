suppressMessages({
    library(posDemux)
    library(glue)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1L]]
bc_cutoff <- as.integer(args[[2L]])

input_freq_table <- glue("results/{sample}/{sample}_frequency_table.txt")
freq_table <- read.table(file = input_freq_table, header = TRUE, sep = "\t")

n_reads <- freq_table$frequency |> sum()

plot_title <- glue("{sample}, n_reads={n_reads}") |> ggtitle()

freq_plot <- freq_plot(freq_table, cutoff = bc_cutoff |> bc_to_freq_cutoff(freq_table = freq_table), type = "density", log_scale_x = TRUE, scale_by_reads = TRUE) + plot_title + theme_bw(base_size = 35)
knee_plot <- knee_plot(freq_table, cutoff = bc_cutoff) + plot_title + theme_bw(base_size = 35)

ggsave(
    filename = "results/{sample}/{sample}_ReadsPerBC.pdf" |> glue(), plot = freq_plot, width = 30,
    height = 20, units = "cm"
)
ggsave(
    filename = "results/{sample}/{sample}_kneePlot.pdf" |> glue(), plot = knee_plot, width = 30,
    height = 20, units = "cm"
)
