# Renames the XS tag to XG it the SAM file because it collides with featureCounts,
# but is still needed by count_genes.py
rule remove_xs_tags:
    input:
        "results/{sample}/{sample}_bwa_{lane}.sam",
    output:
        temp("results/{sample}/{sample}_no_XS_{lane}.sam"),
    shell:
        'sed "s/XS:/XG:/" {input} > {output}'

# Convert to BAM
rule sam_to_bam:
    input:
        "results/{sample}/{sample}_no_XS_{lane}.sam",
    output:
        temp("results/{sample}/{sample}_no_XS_{lane}.bam"),
    shell:
        "samtools view -bS --uncompressed {input} -o {output}"

# Sort before feature counting
rule sam_to_bam_sort:
    input:
        "results/{sample}/{sample}_no_XS_{lane}.bam",
    output:
        temp("results/{sample}/{sample}_sorted_{lane}.bam"),
    shell:
        "samtools sort -u {input} -o {output}"


# Index BAM
rule index_bam:
    input:
        "results/{sample}/{sample}_sorted_{lane}.bam",
    output:
        temp("results/{sample}/{sample}_sorted__{lane}.bam.bai"),
    shell:
        "samtools index {input}"


# Run featureCounts
rule feature_counts:
    input:
        bam="results/{sample}/{sample}_sorted_{lane}.bam",
        bai="results/{sample}/{sample}_sorted_{lane}.bam.bai",
        gff=lambda wildcards: processed_config[wildcards.sample]["annotation"],
    output:
        counts=temp("results/{sample}/{sample}_{lane}.featureCounts.txt"),
        bam=temp("results/{sample}/{sample}_sorted_{lane}.bam.featureCounts.bam"),
    log:
        report="logs/{sample}/featureCounts_{lane}.log",
        # This is really a log file
        counts_summary="logs/{sample}/{sample}_{lane}.featureCounts.txt.summary",
    params:
        feature_tag=lambda wildcards: processed_config[wildcards.sample]["feature_tag"],
        gene_id_attribute=lambda wildcards: processed_config[wildcards.sample]["gene_id_attribute"]
    shell:
        """
        featureCounts -t '{params.feature_tag}' -g {params.gene_id_attribute} -s 1 -a {input.gff} -o {output.counts} -R BAM {input.bam} 2> {log.report}
        mv {output.counts}.summary {log.counts_summary}
        """


# Index featureCounts BAM
rule index_fc_bam:
    input:
        "results/{sample}/{sample}_sorted_{lane}.bam.featureCounts.bam",
    output:
        temp("results/{sample}/{sample}_sorted_{lane}.bam.featureCounts.bam.bai"),
    shell:
        "samtools index {input}"

def get_bam_files_for_sample(wildcards):
    sample = wildcards.sample
    lanes = sample_lanes[sample]
    return [f"results/{sample}/{sample}_sorted_{lane}.bam.featureCounts.bam" for lane in lanes]

def get_bai_files_for_sample(wildcards):
    sample = wildcards.sample
    lanes = sample_lanes[sample]
    return [f"results/{sample}/{sample}_sorted_{lane}.bam.featureCounts.bai" for lane in lanes]

rule count_genes:
    input:
        bam=get_bam_files_for_sample,
        bai=get_bai_files_for_sample,
        barcode_table="results/{sample}/{sample}_selected_barcode_table.sqlite",
    output:
        "results/{sample}/{sample}_gene_count_matrix.txt",
        "results/{sample}/{sample}_umi_count_table.txt"
    log: "logs/{sample}/count_genes.log"
    threads: workflow.cores
    params:
        # This is not the streaming chunk size and is therefore optimized differently
        # Only modify this one if you know what you are doing
        chunk_size=10
    shell:
        "python workflow/scripts/count_genes.py 0 {wildcards.sample} {threads} {params.chunk_size} 2> {log}"
