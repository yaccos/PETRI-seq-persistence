suppressMessages({
    library(Biostrings)
    library(purrr)
    library(glue)
    library(dplyr)
    library(chunked)
    library(DBI)
})

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

args <- commandArgs(trailingOnly = TRUE)

sample  <- args[[1L]]

# sample  <- "random20000"

# bc_cutoff <- 7000L

bc_cutoff  <- as.integer(args[[2L]])
chunk_size  <- as.integer(args[[3L]])
bc_names  <- glue("bc{1:3}")
input_freq_table  <- glue("results/{sample}/{sample}_frequency_table.txt")
input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

log_progress("Reading frequency table")
freq_table <- read.table(file = input_freq_table, header = TRUE, row.names = NULL, sep = "\t")
selected_freq_table <- freq_table[seq_len(min(bc_cutoff, nrow(freq_table))),]

output_database_name <- glue("results/{sample}/{sample}_selected_barcode_table.sqlite")
log_progress("Initializing output barcode database stream")
output_db <- dbConnect(RSQLite::SQLite(), output_database_name)
on.exit(dbDisconnect(output_db), add = TRUE)


log_progress("Selecting reads to keep")
read_table_chunkwise(file = input_table_file, header = TRUE, sep = "\t", chunk_size = chunk_size) |>
inner_join(selected_freq_table)  |> 
mutate(celltag = do.call(paste, c(across(all_of(bc_names)), sep = "_")))  |> 
select(read, UMI, celltag) |> 
write_chunkwise(dbplyr::src_dbi(output_db), "selected_barcodes")


dbExecute(output_db, "CREATE INDEX idx ON selected_barcodes(read)") |> invisible()



log_progress("Writing processed frequency table")
write.table(
    x = selected_freq_table, file = "results/{sample}/{sample}_selected_frequency_table.txt" |> glue(),
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
)

log_progress("DONE")

n_read  <- freq_table$frequency |> sum()
n_selected <- selected_freq_table$frequency |> sum()
selection_percentage <- n_selected / n_read*100
glue("Selected {n_selected} of {n_read} reads ({round(selection_percentage, 2L)}%) based on barcode frequency")
