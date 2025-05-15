def get_input(sample, lane, read):
    sample_config = processed_config[sample]['fastq']
    return sample_config[lane][0 if read == "R1" else 1]

rule fastqc:
    input:
        lambda wildcards: get_input(wildcards.sample, wildcards.lane, wildcards.read),
    output:
        expand(
            "results/{{sample}}/{{sample}}_{{lane}}_{{read}}_fastqc.{format}",
            format= ["html", "zip"]
        ),
    shell:
        # Workaround copy-catted from https://github.com/s-andrews/FastQC/issues/9
        "cutadapt --quiet {input} | fastqc stdin:{wildcards.sample}_{wildcards.lane}_{wildcards.read} -o results/{wildcards.sample}"


rule quality_trim:
    # reverse is reserved for internal use by Snakemake, so we have to rewrite it
    input:
        forward=lambda wildcards: get_input(wildcards.sample, wildcards.lane, "R1"),
        reverse_seq=lambda wildcards: get_input(wildcards.sample, wildcards.lane, "R2"),
    output:
        forward="results/{sample}/{sample}_QF_{lane}_R1.fastq",
        reverse_seq="results/{sample}/{sample}_QF_{lane}_R2.fastq",
    shell:
        "cutadapt -q 10,10 --minimum-length 55:14 --max-n 3 --pair-filter=any -o {output.forward} -p {output.reverse_seq} {input.forward} {input.reverse_seq}"
