rule bwa_index:
    input:
        "{genome}.fa",
    output:
        "{genome}.fa.bwt",
    shell:
        "bwa index {wildcards.genome}.fa"

# Configuration
configfile: "config.yaml"  # Create this file to specify samples, fasta, gff

# Get script directory
SCRIPT_DIR = "PETRI_seq_scripts_v2/scripts/"

# Define final output files
rule all:
    input:
        expand("results/{sample}_v11_threshold_0_matrix.txt", sample=config["samples"])

# BWA alignment
rule bwa_align:
    input:
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=config["fasta"]
    output:
        sai="results/{sample}/{sample}_bwa.sai"
    shell:
        "bwa aln -n 0.06 {input.fasta} {input.fastq} > {output.sai}"

# Convert SAI to SAM
rule bwa_samse:
    input:
        sai="results/{sample}/{sample}_bwa.sai",
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=config["fasta"]
    output:
        sam="results/{sample}/{sample}_bwa.sam"
    shell:
        "bwa samse -n 14 {input.fasta} {input.sai} {input.fastq} > {output.sam}"

# Remove XT tags
rule remove_xt_tags:
    input:
        sam="results/{sample}/{sample}_bwa.sam"
    output:
        sam="results/{sample}/{sample}_no_XT.sam"
    shell:
        "sed \"s/XT:/XN:/\" {input.sam} > {output.sam}"

# Convert to BAM and sort
rule sam_to_bam_sort:
    input:
        sam="results/{sample}/{sample}_no_XT.sam"
    output:
        bam="results/{sample}/{sample}_sorted.bam"
    shell:
        "samtools view -bS {input.sam} | samtools sort - > {output.bam}"

# Index BAM
rule index_bam:
    input:
        bam="results/{sample}/{sample}_sorted.bam"
    output:
        bai="results/{sample}/{sample}_sorted.bam.bai"
    shell:
        "samtools index {input.bam}"

# Run featureCounts
rule feature_counts:
    input:
        bam="results/{sample}/{sample}_sorted.bam",
        bai="results/{sample}/{sample}_sorted.bam.bai",
        gff=config["gff"]
    output:
        counts="results/{sample}_FC/{sample}",
        bam="results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam"
    params:
        outdir=lambda wildcards: f"results/{wildcards.sample}_FC"
    shell:
        """
        mkdir -p {params.outdir}
        featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a {input.gff} -o {output.counts} -R BAM {input.bam}
        """

# Index featureCounts BAM
rule index_fc_bam:
    input:
        bam="results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam"
    output:
        bai="results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam.bai"
    shell:
        "samtools index {input.bam}"

# Add cell barcode
rule add_cell_barcode:
    input:
        bam="results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam",
        bai="results/{sample}_FC/{sample}_sorted.bam.featureCounts.bam.bai"
    output:
        bam="results/{sample}_FC/{sample}_sorted.bam.featureCounts_with_celltag.bam"
    shell:
        "Rscript {SCRIPT_DIR}/add_cell_barcode.R {wildcards.sample}"

# UMI-tools grouping
rule umi_tools_group:
    input:
        bam="results/{sample}_FC/{sample}_sorted.bam.featureCounts_with_celltag.bam"
    output:
        tsv="results/{sample}_FC_directional_grouped_2/{sample}_UMI_counts.tsv",
        bam="results/{sample}_FC_directional_grouped_2/{sample}_group_FC.bam"
    params:
        outdir=lambda wildcards: f"results/{wildcards.sample}_FC_directional_grouped_2"
    shell:
        """
        mkdir -p {params.outdir}
        umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=BX \
        -I {input.bam} \
        --group-out={output.tsv} \
        --method=directional --output-bam -S {output.bam}
        """

# Convert BAM to SAM
rule bam_to_sam:
    input:
        bam="results/{sample}_FC_directional_grouped_2/{sample}_group_FC.bam"
    output:
        sam="results/{sample}_FC_directional_grouped_2/{sample}_group_FC.sam"
    shell:
        "samtools view {input.bam} > {output.sam}"

# Process SAM file
rule process_sam:
    input:
        sam="results/{sample}_FC_directional_grouped_2/{sample}_group_FC.sam"
    output:
        umi="results/{sample}_v11_threshold_0_filtered_mapped_UMIs.txt"
    shell:
        "python {SCRIPT_DIR}/sc_sam_processor_11_generic.py 0 {wildcards.sample}"

# Make matrix
rule make_matrix:
    input:
        umi="results/{sample}_v11_threshold_0_filtered_mapped_UMIs.txt"
    output:
        matrix="results/{sample}_v11_threshold_0_matrix.txt"
    shell:
        "python {SCRIPT_DIR}/make_matrix_mixed_species.py {wildcards.sample}_v11_threshold_0"
