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
    threads: workflow.cores
    shell:
        "pear -f {input.forward} -r {input.reverse_seq} -j {threads} "\
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
    log:
        "logs/{sample}/split_R1_{lane}.log"
    shell:
        "cutadapt -l 58 -o {output} {input} > {log}"


rule split_R2:
    input:
        "results/{sample}/{sample}_QF_{lane}_paired_min75.fastq",
    output:
        temp("results/{sample}/{sample}_QF_{lane}_preR2_paired.fastq"),
    log:
        "logs/{sample}/split_R2_{lane}.log"
    shell:
        "cutadapt -u 58 -o {output} {input} > {log}"


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

# The logic here is that once R1 sequences are used for demultiplexing, they are not longer needed and can be discarded
# before breaking out of the pipeline to determine the BC cutoff
# R2 sequences on the other hand are used past the determination of the BC cutoff and must therefore be perserved for this purpose

def get_lane_files_for_merging(prefix, sample, read):
    template = "{prefix}/{sample}_QF_merged_{lane}_{read}.fastq"
    sample_config = processed_config[sample]
    lanes = sample_config["fastq"].keys()
    return expand(template, lane = lanes, sample= sample, prefix=prefix, read= read)

rule merge_R1:
    input:
        lambda wildcards: get_lane_files_for_merging(wildcards.prefix, wildcards.sample, "R1"),
    output:
        temp("{prefix}/{sample}_QF_merged_R1_all_lanes.fastq"),
    shell:
        "cat {input} > {output}"

rule merge_R2:
    input:
        lambda wildcards: get_lane_files_for_merging(wildcards.prefix, wildcards.sample, "R2"),
    output:
        "{prefix}/{sample}_QF_merged_R2_all_lanes.fastq",
    shell:
        "cat {input} > {output}"
