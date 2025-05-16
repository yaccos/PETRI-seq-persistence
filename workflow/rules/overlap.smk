# Use pear to match read 1 and read 2; for those that overlap, remove reads less than 75bp
rule pear_merge:
    input:
        forward="results/{sample}/{sample}_QF_{lane}_R1.fastq",
        reverse_seq="results/{sample}/{sample}_QF_{lane}_R2.fastq",
    output:
        assembled_reads=temp("results/{sample}/{sample}_QF_{lane}_p.assembled.fastq"),
        unassembled_reads=temp(expand(
            "results/{{sample}}/{{sample}}_QF_{{lane}}_p.unassembled.{direction}.fastq",
            direction=("forward", "reverse"),
        )),
        discarded_reads=temp("results/{sample}/{sample}_QF_{lane}_p.discarded.fastq"),
    log: "logs/{sample}/pear_merge_{lane}.log"
    shell:
        "pear -f {input.forward} -r {input.reverse_seq} "\
        "-o results/{wildcards.sample}/{wildcards.sample}_QF_{wildcards.lane}_p -v 8 -p 0.001 -n 0 > {log}"
        

rule remove_short_reads:
    input:
        "results/{sample}/{sample}_QF_{lane}_p.assembled.fastq",
    output:
        temp("results/{sample}/{sample}_QF_{lane}_paired_min75.fastq"),
    log:
        "logs/{sample}/overlap_removal_{lane}.log"
    shell:
        "cutadapt -m 75 -o {output} {input} > {log}"


rule split_R1:
    input:
        "results/{sample}/{sample}_QF_{lane}_paired_min75.fastq",
    output:
        temp("results/{sample}/{sample}_QF_{lane}_R1_paired.fastq"),
    shell:
        "cutadapt -l 58 -o {output} {input}"


rule split_R2:
    input:
        "results/{sample}/{sample}_QF_{lane}_paired_min75.fastq",
    output:
        temp("results/{sample}/{sample}_QF_{lane}_preR2_paired.fastq"),
    shell:
        "cutadapt -u 58 -o {output} {input}"


rule reverse_compliment_R2:
    input:
        "results/{sample}/{sample}_QF_{lane}_preR2_paired.fastq",
    output:
        temp("results/{sample}/{sample}_QF_{lane}_R2_paired.fastq"),
    shell:
        "seqkit seq -r -p -t DNA {input} -o {output}"


def unassembled_read(wildcards):
    return f"results/{wildcards.sample}/{wildcards.sample}_QF_{wildcards.lane}_p.unassembled.{"forward" if wildcards.read== "R1" else "reverse"}.fastq"


rule merge_reads:
    input:
        assembled="results/{sample}/{sample}_QF_{lane}_{read}_paired.fastq",
        unassembled=unassembled_read,
    output:
        temp("results/{sample}/{sample}_QF_merged_{lane}_{read}.fastq"),
    shell:
        "cat {input.assembled} {input.unassembled} > {output}"

def get_lane_files_for_merging(wildcards):
    template = "{prefix}/{sample}_QF_merged_{lane}_{read}.fastq"
    sample_config = processed_config[wildcards.sample]
    lanes = sample_config["fastq"].keys()
    return expand(template, lane = lanes, sample= wildcards.sample, prefix= wildcards.prefix, read= wildcards.read)

rule merge_lanes:
    input:
        get_lane_files_for_merging,
    output:
        "{prefix}/{sample}_QF_merged_{read}_all_lanes.fastq",
    shell:
        "cat {input} > {output}"
