rule demultiplex:
    input:
        "results/{sample}/{sample}_QF_merged_R1_all_lanes.fastq",
        f"{barcode_dir}/BC1_5p_anchor_v2.fa",
        f"{barcode_dir}/BC2_anchored.fa",
        f"{barcode_dir}/BC3_anchored.fa",
    output:
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
        bc_frame="results/{sample}/{sample}_bc_frame.rds",
        frequency_table="results/{sample}/{sample}_frequency_table.txt",
    shell:
        "Rscript {script_dir}/demultiplexer.R {wildcards.sample}"


rule create_bc_plots:
    input:
        "results/{sample}/{sample}_frequency_table.txt",
    output:
        histogram="results/{sample}/{sample}_ReadsPerBC.pdf",
        knee_plot="results/{sample}/{sample}_kneePlot.pdf",
    params:
        bc_cutoff = lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"]
    shell:
        "Rscript {script_dir}/create_bc_plots.R {wildcards.sample} {params.bc_cutoff}"


rule filter_and_trim:
    input:
        "results/{sample}/{sample}_barcode_table.txt",
        "results/{sample}/{sample}_QF_merged_R2_all_lanes.fastq",
        "results/{sample}/{sample}_bc_frame.rds",
        "results/{sample}/{sample}_frequency_table.txt",
    output:
        trimmed_sequences="results/{sample}/{sample}_2trim.fastq",
        frequency_table="results/{sample}/{sample}_selected_frequency_table.txt",
    params:
        bc_cutoff = lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"]
    shell:
        "Rscript {script_dir}/filter_by_frequency.R {wildcards.sample} {params.bc_cutoff}"
