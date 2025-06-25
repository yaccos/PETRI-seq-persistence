rule diamond_index:
    input:
        "data/{database}.gz",
    output:
        "data/{database}.dmnd",
    params:
        threads=workflow.cores
    log:
        "logs/{database}_diamond_index.log"
    shell:
        "diamond makedb --in {input} -d {output} --threads {params.threads} 2> {log}"

rule diamond:
    input:
        reads="results/{sample}/{sample}_2trim.fastq",
        database=lambda wildcards: processed_config[wildcards.sample]["diamond_database"]
    output:
        "results/{sample}/{sample}_diamond.daa",
    log:
        "logs/{sample}/diamond.log"
    params:
        threads=workflow.cores
    shell:
        "diamond blastx -d {input.database} -q {input.reads} -o {output} --threads {params.threads} -f 100 2> {log}"

rule meganize:
    input:
        daa="results/{sample}/{sample}_diamond.daa",
        database=lambda wildcards: processed_config[wildcards.sample]["megan_database"]
    output:
        "results/{sample}/{sample}_diamond.meganized.daa",
    log:
        "logs/{sample}/meganize.log"
    params:
        threads=workflow.cores
    shell:
        """
        cp {input.daa} {output}
        daa-meganizer -i {output} -mdb {input.database} --threads {params.threads} 2> {log}
        """
