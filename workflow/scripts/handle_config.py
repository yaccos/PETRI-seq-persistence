# The purpose of this script is to untangle the
# Snakemake configuration file and validate its contents

def process_sample(sample_contents, settings, sample_name):
    settings["prefix"] += sample_contents.get("prefix","")
    settings["suffix"] = sample_contents.get("suffix","")
    if "reference_genome" in sample_contents:
        settings["genome"] = settings["prefix"] + sample_contents["reference_genome"] + settings["suffix"]
    elif "genome" not in settings:
        raise ValueError(f"No reference genome specified for sample {sample_name}")
    
    if "reference_annotation" in sample_contents:
        settings["annotation"] = settings["prefix"] + sample_contents["reference_annotation"] + settings["suffix"]
    elif "annotation" not in settings:
        raise ValueError(f"No reference annotation specified for sample {sample_name}")
    
    if "bc_cutoff" in sample_contents:
        # We don't raise an error if this is not defined anywhere
        #  because it is only needed for the latter parts of the workflow
        settings["bc_cutoff"] = sample_contents["bc_cutoff"]
    
    if "feature_tag" in sample_contents:
        settings["feature_tag"] = sample_contents["feature_tag"]
    elif "feature_tag" not in settings:
        raise ValueError(f"No feature tag is defined for sample {sample_name}")
    
    if "gene_id_attribute" in sample_contents:
        settings["gene_id_attribute"] = sample_contents["gene_id_attribute"]
    elif "gene_id_attribute" not in settings:
        raise ValueError(f"No gene identifier attribute is defined for sample {sample_name}")
    
    if "streaming_chunk_size" in sample_contents:
        settings["chunk_size"] = sample_contents["streaming_chunk_size"]
    elif "chunk_size" not in settings:
        raise ValueError(f"No chunk size for FASTQ streaming is defined for sample {sample_name}")
    
    settings["forward_prefix"] = settings["prefix"] + sample_contents.get("forward_prefix", "")
    settings["reverse_prefix"] = settings["prefix"] + sample_contents.get("reverse_prefix", "")

    settings["forward_suffix"] = sample_contents.get("forward_suffix", "") + settings["suffix"]
    settings["reverse_suffix"] = sample_contents.get("reverse_suffix", "") + settings["suffix"]

    lanes = sample_contents.get("lanes", {"": ""})
    if not isinstance(lanes, dict):
        raise ValueError("The lanes field must be a key-value of lane names and identifiers in the input files")
    for lane in lanes:
        if '_' in lane:
            raise ValueError("Lane names are not allowed to have underscores")
        
    settings["fastq"] = {name: (settings["forward_prefix"] + ID + settings["forward_suffix"],
                           settings["reverse_prefix"] + ID + settings["reverse_suffix"]) for name, ID in lanes.items()}
    
    return settings


def config_transform(config):
    settings = {}
    settings["prefix"] = config.get("prefix","")

    if "streaming_chunk_size" in config:
        settings["chunk_size"] = config["streaming_chunk_size"]
    if "reference_genome" in config:
        settings["genome"] = settings["prefix"] + config["reference_genome"]
    if "reference_annotation" in config:
        settings["annotation"] = settings["prefix"] + config["reference_annotation"]

    for key in ["bc_cutoff", "feature_tag", "gene_id_attribute"]:
        if key in config:
            settings[key] = config[key]
    if "samples" not in config:
        raise ValueError("The config file must specify the samples")
    if not isinstance(config["samples"], dict):
        raise ValueError("The sample names must be the keys of the sample listing")
    sample_data = {sample_name: process_sample(contents, settings.copy(), sample_name) for sample_name, contents in config["samples"].items()}
    return sample_data
