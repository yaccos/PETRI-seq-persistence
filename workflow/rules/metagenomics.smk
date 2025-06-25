rule diamond_index:
    input:
        "data/{database}.gz",
    output:
        "data/{database}.dmnd",
    threads: 5
    log:
        "logs/{database}_diamond_index.log"
    shell:
        "diamond makedb --in {input} -d {output} --threads {threads} 2> {log}"

rule diamond:
    input:
        reads="results/{sample}/{sample}_2trim.fastq",
        database=lambda wildcards: processed_config[wildcards.sample]["diamond_database"]
    output:
        "results/{sample}/{sample}_diamond.daa",
    log:
        "logs/{sample}/diamond.log"
    threads: 5
    shell:
        "diamond blastx -d {input.database} -q {input.reads} -o {output} --threads {threads} -f 100 2> {log}"

rule meganize:
    input:
        daa="results/{sample}/{sample}_diamond.daa",
        database=lambda wildcards: processed_config[wildcards.sample]["megan_database"]
    output:
        "results/{sample}/{sample}_diamond.meganized.daa",
    log:
        "logs/{sample}/meganize.log"
    threads: 5
    shell:
        """
        cp {input.daa} {output}
        daa-meganizer -i {output} -mdb {input.database} --threads {threads} 2> {log}
        """
