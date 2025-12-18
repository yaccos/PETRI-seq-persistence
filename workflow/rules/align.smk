from pathlib import Path

bwa_index_extensions = [".bwt", ".amb", ".pac", ".sa", ".ann"]

genome_dict = {}

for sample in sample_names:
    genome_path = processed_config[sample]["genome"]
    genome_name = Path(genome_path).stem
    genome_dict[genome_path] = genome_name


rule bwa_index:
    input:
        "{genome}",
    output:
        temp(multiext("{genome}", *bwa_index_extensions)),
    log:
        lambda wildcards: f"logs/{genome_dict[wildcards.genome]}_index.log"
    conda:
        "../envs/bwa.yml"
    shell:
        "bwa index {input} 2> {log}"

get_bwa_index = lambda wildcards: multiext(processed_config[wildcards.sample]["genome"], *bwa_index_extensions)

# BWA alignment
rule bwa_align:
    input: 
        reads="results/{sample}/{sample}_QF_{lane}_R2.fastq",
        genome=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index
    output:
        temp("results/{sample}/{sample}_bwa_{lane}.sam"),
    log:
        "logs/{sample}/bwa_mem_{lane}.log"
    threads:
        # Ensures there is a core left to do demultiplexing
        max(1, workflow.cores - 1)
    conda:
        "../envs/bwa.yml"
    shell:
        "bwa mem -a -k 19 -Y -t {threads} {input.genome}  {input.reads} > {output} 2> {log}"
