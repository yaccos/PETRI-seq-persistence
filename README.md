# Petrisnake: A secondary analysis pipeline for PETRI-seq data.

This is a Snakemake pipeline for the secondary computational analysis of single cell RNA-seq data from the PETRI-seq protocol (https://www.nature.com/articles/s41564-020-0729-6 and https://www.nature.com/articles/s41586-024-08124-2).

![Workflow rulegraph](./images/rulegraph.svg)

# Dependencies

This pipeline is only tested under Linux running on x86-64 architecture, although it might work on macOS and other CPU architectures if some work is put into dependency management. The workflow will complain, but still work if the system dependency `fuse2fs` is not installed. For running the pipeline the following dependencies are required (tested version):

* Snakemake (9.14.5)
* Apptainer (1.4.5)

The configuration tool GUI is a R Shiny application with the following dependencies:

* R (4.5.2)
* shiny (1.12.1)
* shinyBS (0.61.1)
* purrr (1.2.0)
* yaml (2.3.12)

For convenience, the file `environment.yml` describes the conda dependencies for both Snakemake and the configuration GUI. In order to install it, run `conda env create --file environment.yml`. This should create a conda environment named `petrisnake` with the required dependencies.


The `posDemux` R package (https://github.com/yaccos/posDemux) is a keystone dependency of this workflow. This package is scheduled to be released as part of Bioconductor 3.23  For the pipeline itself, the package is available as an image Docker Hub (https://hub.docker.com/repository/docker/yaccos/posdemux/general) and is automatically used by Snakemake in the workflow provided that the argument `--software-deployment-method apptainer` is used.


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

For using the cutoff selection tools, the `posDemux` package must be run outside the workflow. As of writing Bioconductor release 3.23 is not out yet and the package requires a minimum R version of 4.6.0 as per Bioconductor's policies. This makes the package trickier to install, but there are some options:

* Install the devel version of R which is 4.6. You may use the `rig` version manager (https://github.com/r-lib/rig) for this as it allows you to install multiple versions of R and switch between them as needed.
* Download the sources from GitHub and manually edit the minimum version requirements to your installed R version (everything above 4.4 should work).
* Pull the container image with `posDemux` with Docker and use the bundled RStudio Server instance: `docker run -e PASSWORD=<YOUR_PASS> -p 8787:8787 yaccos/posdemux:0.99.8`. Then open `http://localhost:8787/` in your web browser, type `rstudio` as the username and `<YOUR_PASS>` as the password.

For the rules inside the workflow which do not use R scripts, Snakemakes automatically handles the conda environments for these rules such that no manual intervention is required for the dependency management. These conda environment files are located in `workflow/envs`.

# Configuration



# Running the pipeline

When the barcode threshold to use is unknown, we only need to run the pipeline to the end of the demultiplexer. In that case, run
`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> determine_bc_cutoff`.

Then use the posDemux package to run the interactive cutoff selection tool on the file `results/<sample>/<sample>_barcode_table.txt`.

`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> all`

# Details of the pipeline steps

When running `determine_bc_cutoff`, three rules are run, `quality_trim` and `demultiplex`, and `fastqc`. 

## `fastqc`
Completly separately from the other jobs `fastqc` is run on all input files, the results of which are found in `results/<sample>/qc`.

## `quality_trim`
`quality_trim` is a `cutadapt` command filtering out low-quality reads pairs and reads being too short. The rule is invoked once per input file pair and hence gives a pair of output files per lanes. The output filenames are `results/<sample>/<sample>_QF_<lane>_R1.fastq` and `results/<sample>/<sample>_QF_<lane>_R2.fastq`. The `R1` files are fed into the demultiplexer and are typically never seen by the user, whereas the `R2` files are to be alignment and are hence kept until the `all` rules completes the alignment.

## `demultiplex`
In this step the `posDemux` package demultiplexes the forward reads for each sample. If the sample has multiple lanes, the demultiplexer sequentially streams all the lane files in the same rule invokation.