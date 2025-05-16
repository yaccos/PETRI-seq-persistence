bwa_index_extensions = [".bwt", ".amb", ".pac", ".sa", ".ann"]

rule bwa_index:
    input:
        "data/{genome}.fa",
    output:
        temp(multiext("data/{genome}.fa", *bwa_index_extensions)),
    log:
        "logs/{genome}_index.log"
    shell:
        "bwa index {input} 2> {log}"

get_bwa_index = lambda wildcards: multiext(processed_config[wildcards.sample]["genome"], *bwa_index_extensions)

# BWA alignment
rule bwa_align:
    input:
        reads="results/{sample}/{sample}_2trim.fastq",
        genome=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa.sai"),
    shell:
        "bwa aln -n 0.06 {input.genome} {input.reads} > {output}"


# Convert SAI to SAM
rule bwa_samse:
    input:
        sai="results/{sample}/{sample}_bwa.sai",
        reads="results/{sample}/{sample}_2trim.fastq",
        genome=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa.sam"),
    log:
        "logs/{sample}/bwa_samse.log"
    shell:
        "bwa samse -n 14 {input.genome} {input.sai} {input.reads} > {output} 2> {log}"
