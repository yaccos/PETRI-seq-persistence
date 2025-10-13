# Removes the XT tag from the sam file because it reportably interfer with featureCounts
rule remove_xt_tags:
    input:
        "results/{sample}/{sample}_bwa.sam",
    output:
        temp("results/{sample}/{sample}_no_XT.sam"),
    shell:
        'sed "s/XT:/XN:/" {input} > {output}'


# Convert to BAM and sort
rule sam_to_bam_sort:
    input:
        "results/{sample}/{sample}_no_XT.sam",
    output:
        temp("results/{sample}/{sample}_sorted.bam"),
    shell:
        "samtools view -bS {input} | samtools sort - > {output}"


# Index BAM
rule index_bam:
    input:
        "results/{sample}/{sample}_sorted.bam",
    output:
        temp("results/{sample}/{sample}_sorted.bam.bai"),
    shell:
        "samtools index {input}"


# Run featureCounts
rule feature_counts:
    input:
        bam="results/{sample}/{sample}_sorted.bam",
        bai="results/{sample}/{sample}_sorted.bam.bai",
        gff=lambda wildcards: processed_config[wildcards.sample]["annotation"],
    output:
        counts=temp("results/{sample}/{sample}.featureCounts.txt"),
        bam=temp("results/{sample}/{sample}_sorted.bam.featureCounts.bam"),
    log:
        report="logs/{sample}/{sample}_featureCounts.log",
        # This is really a log file
        counts_summary="logs/{sample}/{sample}.featureCounts.txt.summary",
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
        "results/{sample}/{sample}_sorted.bam.featureCounts.bam",
    output:
        temp("results/{sample}/{sample}_sorted.bam.featureCounts.bam.bai"),
    shell:
        "samtools index {input}"

rule count_genes:
    input:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts.bam",
        bai="results/{sample}/{sample}_sorted.bam.featureCounts.bam.bai",
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
    output:
        "results/{sample}/{sample}_gene_count_matrix.txt"
    log: "logs/{sample}/{sample}_count_genes.log"
    threads: workflow.cores
    params:
        chunk_size=10
    shell:
        "python {script_dir}/count_genes.py 0 {wildcards.sample} {threads} {params.chunk_size} 2> {log}"
