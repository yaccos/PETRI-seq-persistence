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
        "cutadapt -q 10,10 --minimum-length 58:14 --max-n 3 --cores={threads} --pair-filter=any -o {output.forward} -p {output.reverse_seq} {input.forward} {input.reverse_seq} > {log}"

def get_lane_files_for_merging(prefix, sample, read):
    template = "{prefix}/{sample}_QF_{lane}_{read}.fastq"
    lanes = sample_lanes[sample]
    return expand(template, lane=lanes, sample=sample, prefix=prefix, read=read)

rule merge_R1:
    input:
        lambda wildcards: get_lane_files_for_merging(wildcards.prefix, wildcards.sample, "R1"),
    output:
        temp("{prefix}/{sample}_QF_R1_all_lanes.fastq"),
    shell:
        "cat {input} > {output}"
