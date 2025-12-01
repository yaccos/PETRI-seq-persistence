rule demultiplex:
    input:
        "results/{sample}/{sample}_QF_R1_all_lanes.fastq",
        f"{barcode_dir}/BC1_5p_anchor_v2.fa",
        f"{barcode_dir}/BC2_anchored.fa",
        f"{barcode_dir}/BC3_anchored.fa",
    output:
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
        bc_frame="results/{sample}/{sample}_bc_frame.rds",
        freq_table="results/{sample}/{sample}_frequency_table.txt"
    params:
        chunk_size = lambda wildcards: processed_config[wildcards.sample]["chunk_size"]
    log:
        "logs/{sample}/demultiplex.log"
    shell:
        "Rscript \"{script_dir}/demultiplexer.R\" {wildcards.sample} {params.chunk_size} &> {log}"

rule create_bc_plots:
    input:
        "results/{sample}/{sample}_frequency_table.txt",
    output:
        histogram="results/{sample}/{sample}_ReadsPerBC.pdf",
        knee_plot="results/{sample}/{sample}_kneePlot.pdf",
    params:
        bc_cutoff = lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"]
    log:
        "logs/{sample}/bc_plots.log"
    shell:
        "Rscript \"{script_dir}/create_bc_plots.R\" {wildcards.sample} {params.bc_cutoff} &> {log}"


rule select_reads:
    input:
        "results/{sample}/{sample}_frequency_table.txt",
        "results/{sample}/{sample}_barcode_table.txt"
    output:
        "results/{sample}/{sample}_selected_frequency_table.txt",
        temp("results/{sample}/{sample}_selected_barcode_table.sqlite")
    log: 
        "logs/{sample}/select_reads.log"
    params:
        bc_cutoff = lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"],
        chunk_size = lambda wildcards: processed_config[wildcards.sample]["chunk_size"]
    shell:
        "Rscript \"{script_dir}/select_reads_to_keep.R\" {wildcards.sample} {params.bc_cutoff} {params.chunk_size} &> {log}"
        