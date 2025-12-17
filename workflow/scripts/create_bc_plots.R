log_file = snakemake@log[[1]]
log_handle  <- file(log_file, open = "w")
sink(log_handle, append = TRUE, type = "output")
sink(log_handle, append = TRUE, type = "message")

suppressMessages({
    library(posDemux)
    library(glue)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

sample <- snakemake@wildcards[["sample"]]
bc_cutoff <- snakemake@params[["bc_cutoff"]]

input_freq_table <- snakemake@input[[1]]
freq_table <- read.table(file = input_freq_table, header = TRUE, sep = "\t")

n_reads <- freq_table$frequency |> sum()

plot_title <- glue("{sample}, n_reads={n_reads}") |> ggtitle()

freq_plot <- freq_plot(freq_table, cutoff = bc_cutoff |> bc_to_freq_cutoff(freq_table = freq_table), type = "density", log_scale_x = TRUE, scale_by_reads = TRUE) + plot_title + theme_bw(base_size = 35)
knee_plot <- knee_plot(freq_table, cutoff = bc_cutoff) + plot_title + theme_bw(base_size = 35)

ggsave(
    filename = snakemake@output[["histogram"]], plot = freq_plot, width = 30,
    height = 20, units = "cm"
)
ggsave(
    filename = snakemake@output[["knee_plot"]], plot = knee_plot, width = 30,
    height = 20, units = "cm"
)

sink(type="message")
sink(type="output")
close(log_handle)
