rule bwa_index:
    input:
        "{genome}.fa",
    output:
        "{genome}.fa.bwt",
    shell:
        "bwa index {wildcards.genome}.fa"
