bwa_index_extensions = [".bwt", ".amb", ".pac", ".sa", ".ann"]

rule bwa_index:
    input:
        "resources/{genome}",
    output:
        temp(multiext("resources/{genome}", *bwa_index_extensions)),
    log:
        "logs/{genome}_index.log"
    shell:
        "bwa index {input} 2> {log}"

get_bwa_index = lambda wildcards: multiext(processed_config[wildcards.sample]["genome"], *bwa_index_extensions)

# BWA alignment
rule bwa_align:
    input: 
        reads="results/{sample}/{sample}_QF_R2_all_lanes.fastq",
        genome=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa.sam"),
    log:
        "logs/{sample}/bwa_mem.log"
    threads:
        # Ensures there is a core left to do demultiplexing
        max(1, workflow.cores - 1)
    shell:
        "bwa mem -a -k 19 -Y -t {threads} {input.genome}  {input.reads} > {output} 2> {log}"
