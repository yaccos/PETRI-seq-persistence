rule fastqc:
    input:
        "data/{sample}/{fastq_name}.fastq.gz",
    output:
        expand(
            "results/{{sample}}/{{fastq_name}}_fastqc.{format}",
            format= ["html", "zip"]
        ),
    shell:
        "fastqc {input} -o results/{sample}"


rule quality_trim:
    # reverse is reserved for internal use by Snakemake, so we have to rewrite it
    input:
        forward="data/{sample}/{sample}_S1_L00{lane}_R1_001.fastq.gz",
        reverse_seq="data/{sample}/{sample}_S1_L00{lane}_R2_001.fastq.gz",
    output:
        forward="results/{sample}/{sample}_QF_L00{lane}_R1_001.fastq",
        reverse_seq="results/{sample}/{sample}_QF_L00{lane}_R2_001.fastq",
    shell:
        f"cutadapt -q 10,10 --minimum-length 55:14 --max-n 3 --pair-filter=any -o {output.forward} -p {output.reverse_seq} {input.forward} {input.reverse_seq}"
