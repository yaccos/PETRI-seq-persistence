rule fastqc:
    input:
        lambda wildcards: get_fastqc_input(wildcards.sample, wildcards.filename),
    output:
        expand(
            "results/{{sample}}/qc/{{filename}}_fastqc.{format}",
            format= ["html", "zip"]
        ),
    log:
        "logs/{sample}/qc/{filename}_fastqc.log"
    shell:
        "fastqc {input} -o results/{wildcards.sample}/qc 2> {log}"


rule quality_trim:
    # reverse is reserved for internal use by Snakemake, so we have to rewrite it
    input:
        forward=lambda wildcards: get_input(wildcards.sample, wildcards.lane, "R1"),
        reverse_seq=lambda wildcards: get_input(wildcards.sample, wildcards.lane, "R2"),
    output:
        forward=temp("results/{sample}/{sample}_QF_{lane}_R1.fastq"),
        reverse_seq=temp("results/{sample}/{sample}_QF_{lane}_R2.fastq"),
    # For this task, cutadapt appears not to take advantage of more than 3 cores 
    threads: min(3, workflow.cores)
    log:
        "logs/{sample}/QF_{lane}.log"
    shell:
        "cutadapt -q 10,10 --minimum-length 58:14 --max-n 3 --cores={threads} --pair-filter=any -g NNNNNNNNAGAATACACGACGCTCTTCCGATCT -o {output.forward} -p {output.reverse_seq} {input.forward} {input.reverse_seq} > {log}"

def get_lane_files_for_merging(prefix, sample, read):
    template = "{prefix}/{sample}_QF_{lane}_{read}.fastq"
    sample_config = processed_config[sample]
    lanes = sample_config["fastq"].keys()
    return expand(template, lane=lanes, sample=sample, prefix=prefix, read=read)

rule merge_R1:
    input:
        lambda wildcards: get_lane_files_for_merging(wildcards.prefix, wildcards.sample, "R1"),
    output:
        temp("{prefix}/{sample}_QF_R1_all_lanes.fastq"),
    shell:
        "cat {input} > {output}"


# The logic here is that once R1 sequences are used for demultiplexing, they are not longer needed and can be discarded
# before breaking out of the pipeline to determine the BC cutoff
# R2 sequences on the other hand are used past the determination of the BC cutoff and must therefore be perserved for this purpose

rule merge_R2:
    input:
        lambda wildcards: get_lane_files_for_merging(wildcards.prefix, wildcards.sample, "R2"),
    output:
        "{prefix}/{sample}_QF_R2_all_lanes.fastq",
    shell:
        "cat {input} > {output}"
