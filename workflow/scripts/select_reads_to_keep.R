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
input_frequency_table  <- glue("results/{sample}/{sample}_frequency_table.txt")
input_table_file <- glue("results/{sample}/{sample}_barcode_table.txt")

log_progress("Reading frequency table")
freq_table <- read.table(file = input_frequency_table, header = TRUE, row.names = NULL, sep = "\t")
selected_freq_table <- freq_table[seq_len(min(bc_cutoff, nrow(freq_table))),]
log_progress("Inititializing input barcode table stream")
input_con <- file(input_table_file, open = "r")
on.exit(close(input_con), add=TRUE)

output_database_name <- glue("results/{sample}/{sample}_selected_barcode_table.sqlite")
log_progress("Initializing output barcode database stream")
output_db <- dbConnect(RSQLite::SQLite(), output_database_name)

dbExecute(output_db, "DROP TABLE IF EXISTS selected_barcodes")
dbExecute(
    output_db,
    "CREATE TABLE selected_barcodes (
        read TEXT NOT NULL PRIMARY KEY,
        UMI TEXT NOT NULL,
        celltag TEXT NOT NULL
    ) WITHOUT ROWID"
)
dbExecute(output_db, "CREATE INDEX idx ON selected_barcodes(read)")

# dbWriteTable(output_db, "barcodes")
on.exit(dbDisconnect(output_db), add = TRUE)



log_progress("Selecting reads to keep")
read_table_chunkwise(file = input_table_file, header = TRUE, sep = "\t", chunk_size = chunk_size) |> 
(\(table) filter(table, posDemux::row_match(table, selected_freq_table)))() |> 
mutate(read, UMI, pasted_celltag = rlang::exec(!!! bc_names), .keep= "none") |> 
write_chunkwise()


header <- read.table(input_con, nrows = 0, header=TRUE,
 colClasses = c("character","character","character"))
repeat {
    chunk <- read.table(
        input_con,
        nrows = chunk_size,
        header = FALSE,
        col.names = names(header),
        colClasses = c("character","character","character")
    )
    if (!nrow(chunk)) break

    common_rows <- posDemux::row_match(chunk, selected_freq_table)

    pasted_celltag <- chunk  |>
     filter(common_rows)  |>
     select(!!! bc_names) |>
      (\(x) do.call(paste, c(x, sep = "_")))()

    
    selected_barcode_table  <- data.frame(read=chunk$read[common_rows], UMI=chunk$UMI[common_rows],
    celltag=pasted_celltag)
    dbAppendTable(output_db, "selected_barcodes", selected_barcode_table)
}





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
