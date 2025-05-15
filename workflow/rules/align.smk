bwa_index_extensions = [".bwt", ".amb", ".pac", ".sa", ".ann"]

rule bwa_index:
    input:
        "{genome}.fa",
    output:
        temp(multiext("{genome}.fa", *bwa_index_extensions)),
    shell:
        "bwa index {wildcards.genome}.fa"

get_bwa_index = lambda wildcards: multiext(processed_config[wildcards.sample]["genome"], *bwa_index_extensions)

# BWA alignment
rule bwa_align:
    input:
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa.sai"),
    shell:
        "bwa aln -n 0.06 {input.fasta} {input.fastq} > {output}"


# Convert SAI to SAM
rule bwa_samse:
    input:
        sai="results/{sample}/{sample}_bwa.sai",
        fastq="results/{sample}/{sample}_2trim.fastq",
        fasta=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa.sam"),
    shell:
        "bwa samse -n 14 {input.fasta} {input.sai} {input.fastq} > {output}"
