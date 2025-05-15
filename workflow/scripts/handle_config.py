# The purpose of this script is to untangle the
# Snakemake configuration file and validate its contents

def check_lane(lane):
    if '_' in lane:
        raise ValueError("Lane names are not allowed to have underscores")
    return lane

def process_sample(sample_contents, settings, sample_name):
    settings["prefix"] += sample_contents.get("prefix","")
    settings["suffix"] += sample_contents.get("suffix","")
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
    
    settings["forward_prefix"] = settings["prefix"] + sample_contents.get("forward_prefix", "")
    settings["reverse_prefix"] = settings["prefix"] + sample_contents.get("reverse_prefix", "")

    settings["forward_suffix"] = sample_contents.get("forward_suffix", "") + settings["suffix"]
    settings["reverse_suffix"] = sample_contents.get("reverse_suffix", "") + settings["suffix"]

    lanes = sample_contents.get("lanes", [""])
    if not isinstance(lanes, list):
        raise ValueError("The lanes field must be a simple sequence of lane names")
    for lane in lanes:
        if '_' in lane:
            raise ValueError("Lane names are not allowed to have underscores")
        
    settings["fastq"] = {check_lane(lane): (settings["forward_prefix"] + lane + settings["forward_suffix"],
                           settings["reverse_prefix"] + lane + settings["reverse_suffix"]) for lane in lanes}
    
    return settings


def config_transform(config):
    settings = {}
    settings["prefix"] = config.get("prefix","")
    settings["suffix"] = config.get("suffix", "")
    if "reference_genome" in config:
        settings["genome"] = settings["prefix"] + config["reference_genome"] + settings["suffix"]
    if "reference_annotation" in config:
        settings["annotation"] = settings["prefix"] + config["reference_annotation"] + settings["suffix"]
    if "bc_cutoff" in config:
        settings["bc_cutoff"] = config["bc_cutoff"]
    if "samples" not in config:
        raise ValueError("The config file must specify the samples")
    if not isinstance(config["samples"], dict):
        raise ValueError("The sample names must be the keys of the sample listing")
    sample_data = {sample_name: process_sample(contents, settings.copy(), sample_name) for sample_name, contents in config["samples"].items()}
    return sample_data
