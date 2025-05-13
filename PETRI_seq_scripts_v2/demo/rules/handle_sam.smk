# Removes the XT tag from the sam file because it reportably interfer with featureCount
rule remove_xt_tags:
    input:
        sam="results/{sample}/{sample}_bwa.sam",
    output:
        sam="results/{sample}/{sample}_no_XT.sam",
    shell:
        'sed "s/XT:/XN:/" {input.sam} > {output.sam}'


# Convert to BAM and sort
rule sam_to_bam_sort:
    input:
        sam="results/{sample}/{sample}_no_XT.sam",
    output:
        bam="results/{sample}/{sample}_sorted.bam",
    shell:
        "samtools view -bS {input.sam} | samtools sort - > {output.bam}"


# Index BAM
rule index_bam:
    input:
        bam="results/{sample}/{sample}_sorted.bam",
    output:
        bai="results/{sample}/{sample}_sorted.bam.bai",
    shell:
        "samtools index {input.bam}"


# Run featureCounts
rule feature_counts:
    input:
        bam="results/{sample}/{sample}_sorted.bam",
        bai="results/{sample}/{sample}_sorted.bam.bai",
        gff=reference_annotation,
    output:
        counts="results/{sample}/{sample}.featureCounts.txt",
        bam="results/{sample}/{sample}_sorted.bam.featureCounts.bam",
    params:
        outdir=lambda wildcards: f"results/{wildcards.sample}_FC",
    shell:
        """
        featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a {input.gff} -o {output.counts} -R BAM {input.bam}
        """


# Index featureCounts BAM
rule index_fc_bam:
    input:
        bam="results/{sample}_sorted.bam.featureCounts.bam",
    output:
        bai="results/{sample}_sorted.bam.featureCounts.bam.bai",
    shell:
        "samtools index {input.bam}"


# Add cell barcode
rule add_cell_barcode:
    input:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts.bam",
        bai="results/{sample}/{sample}_sorted.bam.featureCounts.bam.bai",
        barcode_table="results/{sample}_barcode_table.txt",
    output:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam",
    shell:
        "Rscript {script_dir}/add_cell_barcode.R {wildcards.sample}"


# UMI-tools grouping
rule umi_tools_group:
    input:
        bam="results/{sample}/{sample}_sorted.bam.featureCounts_with_celltag.bam",
    output:
        tsv="results/{sample}/{sample}_UMI_counts.tsv",
        bam="results/{sample}/{sample}_group_FC.bam",
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
        bam="results/{sample}/{sample}_group_FC.bam",
    output:
        sam="results/{sample}/{sample}_group_FC.sam",
    shell:
        "samtools view {input.bam} > {output.sam}"


# Process SAM file
rule process_sam:
    input:
        sam="results/{sample}/{sample}_group_FC.sam",
    output:
        umi="results/{sample}_filtered_mapped_UMIs.txt",
    shell:
        "python {script_dir}/sc_sam_processor.py 0 {wildcards.sample}"


# Make matrix
rule make_matrix:
    input:
        umi="results/{sample}_filtered_mapped_UMIs.txt",
    output:
        matrix="results/{sample}_gene_count_matrix.txt",
    shell:
        "python {script_dir}/make_matrix_mixed_species.py {wildcards.sample}"
