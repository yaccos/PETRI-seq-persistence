def get_input(sample, lane, read):
    sample_config = processed_config[sample]['fastq']
    return sample_config[lane][0 if read == "R1" else 1]

def get_path_stem(path: str):
    # In case we get compressed files, they often have double extensions,
    # so we remove the final .gz in that case. Otherwise, only the final extension is removed
    uncompressed_path = path.removesuffix(".gz")
    stem = pathlib.Path(uncompressed_path).stem
    return stem

def expand_fastqc(sample):
    sample_config = processed_config[sample]
    template = "results/{sample}/qc/{filename}_fastqc.{format}"
    input_lanes = sample_config["fastq"].values()
    input_filestems = [get_path_stem(filename) for lane in input_lanes for filename in lane]
    return expand(
    template,
    sample=sample,
    filename=input_filestems,
    format=["html", "zip"],
    )

def get_fastqc_input(sample, filename):
    sample_config = processed_config[sample]
    for lane_files in sample_config["fastq"].values():
        for path in lane_files:
            if get_path_stem(path) == filename:
                return path
    raise ValueError("FastQC input file not found")

def expand_QF_files(sample):
    return [f"results/{sample}/{sample}_QF_{lane}_R2.fastq" for sample in sample_names for lane in sample_lanes[sample]]


def expand_sample(sample, file):
    return expand("results/{sample}/{sample}_{file}", sample=sample_names, file=file)

def intermediate_files():
    samples = processed_config.keys()
    templates = ["results/{sample}/{sample}_barcode_table.txt", "results/{sample}/{sample}_bc_frame.rds",
     "results/{sample}/{sample}_frequency_table.txt"]
    files = [expand(template, sample=samples) for template in templates]
    return files
