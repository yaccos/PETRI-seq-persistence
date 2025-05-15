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
        fasta=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=lambda wildcards: processed_config[wildcards.sample]["genome"] + ".bwt"
    output:
        "results/{sample}/{sample}_bwa.sai",
    shell:
        "bwa aln -n 0.06 {input.fasta} {input.fastq} > {output}"


# Convert SAI to SAM
rule bwa_samse:
    input:
        sai="results/{sample}/{sample}_bwa.sai",
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=lambda wildcards: processed_config[wildcards.sample]["genome"],
    output:
        "results/{sample}/{sample}_bwa.sam",
    shell:
        "bwa samse -n 14 {input.fasta} {input.sai} {input.fastq} > {output}"
