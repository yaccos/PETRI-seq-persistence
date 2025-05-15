# Removes the XT tag from the sam file because it reportably interfer with featureCount
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
        # This is really a log file
        counts_summary="results/{sample}/{sample}.featureCounts.txt.summary",
        bam=temp("results/{sample}/{sample}_sorted.bam.featureCounts.bam"),
    shell:
        """
        featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a {input.gff} -o {output.counts} -R BAM {input.bam}
        """


# Index featureCounts BAM
rule index_fc_bam:
    input:
        "results/{sample}/{sample}_sorted.bam.featureCounts.bam",
    output:
        temp("results/{sample}/{sample}_sorted.bam.featureCounts.bam.bai"),
    shell:
        "samtools index {input}"


# Add cell barcode
rule add_cell_barcode:
    input:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts.bam",
        bai="results/{sample}/{sample}_sorted.bam.featureCounts.bam.bai",
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
    output:
        temp("results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam"),
        temp("results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam.bai")
    shell:
        "Rscript {script_dir}/add_cell_barcode.R {wildcards.sample}"


# UMI-tools grouping
rule umi_tools_group:
    input:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam",
        bai="results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam.bai"
    output:
        tsv=temp("results/{sample}/{sample}_UMI_counts.tsv"),
        bam=temp("results/{sample}/{sample}_group_FC.bam"),
    shell:
        """
        umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=BX \
        -I {input.bam} \
        --group-out={output.tsv} \
        --method=directional --output-bam -S {output.bam}
        """


# Convert BAM to SAM
rule bam_to_sam:
    input:
        "results/{sample}/{sample}_group_FC.bam",
    output:
        temp("results/{sample}/{sample}_group_FC.sam"),
    shell:
        "samtools view {input} > {output}"


# Process SAM file
rule process_sam:
    input:
        "results/{sample}/{sample}_group_FC.sam",
    output:
        temp("results/{sample}/{sample}_filtered_mapped_UMIs.txt"),
    shell:
        "python {script_dir}/sc_sam_processor.py 0 {wildcards.sample}"


# Make matrix
rule make_matrix:
    input:
        "results/{sample}/{sample}_filtered_mapped_UMIs.txt",
    output:
        "results/{sample}/{sample}_gene_count_matrix.txt",
    shell:
        "python {script_dir}/make_matrix_mixed_species.py {wildcards.sample}"
