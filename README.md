# Petrisnake: A secondary analysis pipeline for PETRI-seq data.

This is a Snakemake pipeline for the secondary computational analysis of single cell RNA-seq data from the PETRI-seq protocol (https://www.nature.com/articles/s41564-020-0729-6 and https://www.nature.com/articles/s41586-024-08124-2).

![Workflow rulegraph](./images/rulegraph.svg)

# Dependencies

This pipeline is only tested under Linux running on x86-64 architecture, although it might work on macOS and other CPU architectures. The workflow will complain, but still work if the system dependency `fuse2fs` is not installed. The `posDemux` R package (https://github.com/yaccos/posDemux) is a keystone dependency of this workflow. This package is scheduled to be released as part of  For the pipeline itself, the package is available as an image on 

 For using the cutoff selection tools,

When running the container, you may get the warning messages:
```
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
6: Setting LC_PAPER failed, using "C" 
7: Setting LC_MEASUREMENT failed, using "C"
```

These warnings can be safely ignored, but if you still want to silence them, run `export LC_ALL=C` in your shell before running Snakemake.

# Configuration

# Running the pipeline

When the barcode threshold to use is unknown, we only need to run the pipeline to the end of the demultiplexer. In that case, run
`snakemake --software-deployment-method apptainer --cores <NUMBER_OF_CORES> determine_bc_cutoff`.

Then use the posDemux package to run the interactive cutoff selection tool on the file `results/<sample>/<sample>_barcode_table.txt`.

`snakemake --software-deployment-method apptainer --cores <NUMBER_OF_CORES> all`
