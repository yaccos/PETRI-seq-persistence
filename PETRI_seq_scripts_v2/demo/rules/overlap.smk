# Use pear to match read 1 and read 2; for those that overlap, remove reads less than 75bp


rule pear_merge:
    input:
        forward="results/{sample}/{sample}_QF_L00{lane}_R1_001.fastq",
        reverse_seq="results/{sample}/{sample}_QF_L00{lane}_R2_001.fastq",
    output:
        assembled_reads="results/{sample}/{sample}_QF_L00{lane}_p.assembled.fastq",
        unassembled_reads=expand(
            "results/{{sample}}/{{sample}}_QF_L00{{lane}}_p.unassembled.{direction}.fastq",
            direction=("forward", "reverse"),
        ),
        discarded_reads="results/{sample}/{sample}_QF_L00{lane}_p.discarded.fastq",
    shell:
        "pear -f {input.forward} -r {input.reverse_seq} -o results/{wildcards.sample}/{wildcards.sample}_QF_L00{wildcards.lane}_p -v 8 -p 0.001 -n 0"


rule remove_short_reads:
    input:
        "results/{sample}/{sample}_QF_L00{lane}_p.assembled.fastq",
    output:
        "results/{sample}/{sample}_QF_L00{lane}_paired_min75_001.fastq",
    shell:
        "cutadapt -m 75 -o {output} {input}"


rule split_R1:
    input:
        "results/{sample}/{sample}_QF_L00{lane}_paired_min75_001.fastq",
    output:
        "results/{sample}/{sample}_QF_L00{lane}_R1_paired.fastq",
    shell:
        "cutadapt -l 58 -o {output} {input}"


rule split_R2:
    input:
        "results/{sample}/{sample}_QF_L00{lane}_paired_min75_001.fastq",
    output:
        "results/{sample}/{sample}_QF_L00{lane}_preR2_paired.fastq",
    shell:
        "cutadapt -u 58 -o {output} {input}"


rule reverse_compliment_R2:
    input:
        "results/{sample}/{sample}_QF_L00{lane}_preR2_paired.fastq",
    output:
        "results/{sample}/{sample}_QF_L00{lane}_R2_paired.fastq",
    shell:
        "seqkit seq -r -p -t DNA {input} -o {output}"


def unassembled_read(wildcards):
    return f"results/{wildcards.sample}/{wildcards.sample}_QF_L00{wildcards.lane}_p.unassembled.{"forward" if wildcards.read== "R1" else "reverse"}.fastq"


rule merge_reads:
    input:
        assembled="results/{sample}/{sample}_QF_L00{lane}_{read}_paired.fastq",
        unassembled=unassembled_read,
    output:
        "results/{sample}/{sample}_QF_merged_L00{lane}_{read}.fastq",
    shell:
        "cat {input.assembled} {input.unassembled} > {output}"


rule merge_lanes:
    input:
        [f"{{prefix}}_L00{lane}_{{read}}.fastq" for lane in range(1, n_lanes + 1)],
    output:
        "{prefix}_{read}_all_lanes.fastq",
    shell:
        "cat {input} > {output}"
