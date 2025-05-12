bc_cutoff = 7000


rule demultiplex:
    input:
        f"results/{sample}/{sample}_QF_merged_R1_all_lanes.fastq",
        f"{barcode_dir}/BC1_5p_anchor_v2.fa",
        f"{barcode_dir}/BC2_anchored.fa",
        f"{barcode_dir}/BC3_anchored.fa",
    output:
        barcode_table = f"results/{sample}_barcode_table.txt",
        bc_frame=f"results/{sample}_bc_frame.rds",
        frequency_table=f"results/{sample}_frequency_table.txt",
    shell:
        f"Rscript {script_dir}/demultiplexer.R {sample}"


rule create_bc_plots:
    input:
        f"results/{sample}_frequency_table.txt",
    output:
        histogram=f"results/{sample}_ReadsPerBC.pdf",
        knee_plot=f"results/{sample}_kneePlot.pdf",
    shell:
        f"Rscript {script_dir}/create_bc_plots.R {sample} {bc_cutoff}"


rule filter_and_trim:
    input:
        f"results/{sample}_barcode_table.txt",
        f"results/{sample}/{sample}_QF_merged_R2_all_lanes.fastq",
        f"results/{sample}_bc_frame.rds",
        f"results/{sample}_frequency_table.txt",
    output:
        trimmed_sequences=f"results/{sample}/{sample}_2trim.fastq",
        frequency_table=f"results/{sample}_selected_frequency_table.txt",
    shell:
        f"Rscript {script_dir}/filter_by_frequency.R {sample} {bc_cutoff}"
