rule bwa_index:
    input:
        "{genome}.fa",
    output:
        "{genome}.fa.bwt",
    shell:
        "bwa index {wildcards.genome}.fa"


# BWA alignment
rule bwa_align:
    input:
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=reference_genome,
        index=f"{reference_genome}.bwt"
    output:
        sai="results/{sample}/{sample}_bwa.sai",
    shell:
        "bwa aln -n 0.06 {input.fasta} {input.fastq} > {output.sai}"


# Convert SAI to SAM
rule bwa_samse:
    input:
        sai="results/{sample}/{sample}_bwa.sai",
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=reference_genome,
    output:
        sam="results/{sample}/{sample}_bwa.sam",
    shell:
        "bwa samse -n 14 {input.fasta} {input.sai} {input.fastq} > {output.sam}"
